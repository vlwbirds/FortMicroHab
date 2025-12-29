# ============================================================
# 0) Libraries + terra settings
# ============================================================
library(here)
library(terra)
library(sf)

terraOptions(
  todisk = TRUE,
  memfrac = 0.6,
  tempdir = here("terra_tmp")  # create once
)
dir.create(here("terra_tmp"), showWarnings = FALSE, recursive = TRUE)

# Toggle masking later (OFF for now)
USE_SHADOW_MASK <- FALSE

# ============================================================
# 1) Inputs
# ============================================================
aoi_path <- here("data", "fort_boundary.kml")
aoi_v <- terra::vect(aoi_path)
# Project AOI to NAIP CRS later after NAIP is opened

naip_dir <- here("data", "naip_tiles_2023")
dem_dir  <- here("data", "dem")

out_dir_tiles <- here("output", "classified_tiles")
dir.create(out_dir_tiles, showWarnings = FALSE, recursive = TRUE)

out_vrt_class <- here("output", "classified_mosaic.vrt")
out_tif_class <- here("output", "classified_mosaic.tif")

# ============================================================
# 2) Build VRTs (fast, no mosaic yet)
# ============================================================
naip_files <- list.files(naip_dir, pattern="\\.tif(f)?$", full.names=TRUE, ignore.case=TRUE)
stopifnot(length(naip_files) > 0)

dem_files  <- list.files(dem_dir,  pattern="\\.tif(f)?$", full.names=TRUE, ignore.case=TRUE)
stopifnot(length(dem_files) > 0)

naip_vrt <- here("output", "naip.vrt")
dem_vrt  <- here("output", "dem.vrt")
dir.create(here("output"), showWarnings = FALSE, recursive = TRUE)

# Build VRTs (overwrites)
terra::vrt(naip_files, filename = naip_vrt, overwrite = TRUE)
terra::vrt(dem_files,  filename = dem_vrt,  overwrite = TRUE)

naip <- rast(naip_vrt)   # 4-band NAIP mosaic-as-VRT
dem  <- rast(dem_vrt)    # DEM mosaic-as-VRT

# Match AOI CRS to NAIP
aoi_v <- project(aoi_v, crs(naip))

# Optional: crop both to AOI bbox early to reduce read overhead
bb <- ext(aoi_v)
naip <- crop(naip, bb)
dem  <- crop(dem,  bb)

# ============================================================
# 3) Band layout helper: NAIP 4-band -> R,G,B,NIR
# ============================================================
guess_layout <- function(x4) {
  if (nlyr(x4) != 4) stop("Expected 4-band NAIP (RGBN). Got ", nlyr(x4))
  m <- sapply(1:4, function(i) global(x4[[i]], "mean", na.rm=TRUE)[[1]])
  nir_idx <- which.max(m)
  if (nir_idx == 4) return(list(R=1,G=2,B=3,NIR=4))
  if (nir_idx == 1) return(list(NIR=1,R=2,G=3,B=4))
  list(R=1,G=2,B=3,NIR=4)
}

idx <- guess_layout(naip)
naip4 <- c(naip[[idx$R]], naip[[idx$G]], naip[[idx$B]], naip[[idx$NIR]])
names(naip4) <- c("R","G","B","NIR")

# scale to 0-1 if needed
rng <- global(naip4[[1]], c("min","max"), na.rm=TRUE)[1,]
if (rng[2] > 1.5) naip4 <- naip4 / 255

# ============================================================
# 4) Align DEM to NAIP grid (do this ONCE, via VRT resample)
# ============================================================
# If DEM is 1m and NAIP is 0.3m, aligning DEM to NAIP grid is heavy but doable.
# For speed you can instead downsample NAIP to 1m. If you want 30cm classification, keep as is.

dem_rs <- resample(dem, naip4[[1]], method="bilinear")  # DEM now on NAIP grid
names(dem_rs) <- "DEM"

# ============================================================
# 5) Build predictors (global)
# ============================================================
ndvi <- (naip4$NIR - naip4$R) / (naip4$NIR + naip4$R)
ndwi <- (naip4$G   - naip4$NIR) / (naip4$G   + naip4$NIR)
names(ndvi) <- "NDVI"; names(ndwi) <- "NDWI"

# texture proxy: NDVI SD (5x5)
ndvi_sd <- focal(
  ndvi, w=matrix(1,5,5),
  fun=function(x, ...) sd(x, na.rm=TRUE)
)
names(ndvi_sd) <- "NDVI_sd"

# slope/aspect from DEM on NAIP grid
slope_deg  <- terrain(dem_rs, v="slope", unit="degrees")
aspect_rad <- terrain(dem_rs, v="aspect", unit="radians")

thresh <- 0.5
asp_sin <- ifel(is.na(slope_deg) | is.na(aspect_rad) | slope_deg < thresh, 0, sin(aspect_rad))
asp_cos <- ifel(is.na(slope_deg) | is.na(aspect_rad) | slope_deg < thresh, 0, cos(aspect_rad))
names(slope_deg) <- "slope"
names(asp_sin)   <- "asp_sin"
names(asp_cos)   <- "asp_cos"

tpi <- terrain(dem_rs, v="TPI"); names(tpi) <- "TPI"
tri <- terrain(dem_rs, v="TRI"); names(tri) <- "TRI"

preds <- c(naip4, ndvi, ndwi, ndvi_sd, dem_rs, slope_deg, asp_sin, asp_cos, tpi, tri)

# Clip predictors to AOI (mask) to avoid sampling outside
preds <- mask(preds, aoi_v)

# ============================================================
# 6) Fit k-means ONCE (global sample across AOI)
# ============================================================
set.seed(42)
samp <- spatSample(preds, size=150000, method="random", na.rm=TRUE, as.points=FALSE)
X <- na.omit(as.data.frame(samp, cells=FALSE))
k <- 12  # set your endgame here; or run multiple ks if you want

X_scaled <- scale(X)
mu <- attr(X_scaled, "scaled:center")
sdv <- attr(X_scaled, "scaled:scale")
km <- kmeans(X_scaled, centers=k, nstart=10)
centers <- km$centers

assign_cluster <- function(v, centers, mu, sd) {
  if (any(is.na(v))) return(NA_integer_)
  z <- (v - mu) / sd
  d2 <- rowSums((centers - matrix(z, nrow(centers), ncol(centers), byrow=TRUE))^2)
  which.min(d2)
}

# ============================================================
# 7) Tile grid over AOI (in pixels) and run chunked classification
# ============================================================
# Choose chunk size in pixels (NOT 100000 unless you mean it)
tile_ncol <- 10000
tile_nrow <- 10000

# Build tile extents aligned to raster grid
r_template <- preds[[1]]  # a single layer for geometry
e_aoi <- ext(aoi_v)

# helper to snap to raster grid
snap_extent <- function(e, r) {
  resx <- res(r)[1]; resy <- res(r)[2]
  ex <- e
  ex$xmin <- floor((ex$xmin - ext(r)$xmin) / resx) * resx + ext(r)$xmin
  ex$xmax <- ceiling((ex$xmax - ext(r)$xmin) / resx) * resx + ext(r)$xmin
  ex$ymin <- floor((ex$ymin - ext(r)$ymin) / resy) * resy + ext(r)$ymin
  ex$ymax <- ceiling((ex$ymax - ext(r)$ymin) / resy) * resy + ext(r)$ymin
  ex
}

e_aoi <- snap_extent(e_aoi, r_template)

# generate a grid of extents
xs <- seq(e_aoi$xmin, e_aoi$xmax, by = tile_ncol * res(r_template)[1])
ys <- seq(e_aoi$ymin, e_aoi$ymax, by = tile_nrow * res(r_template)[2])

tile_exts <- list()
tid <- 0L
for (x0 in xs) {
  for (y0 in ys) {
    tid <- tid + 1L
    tile_exts[[tid]] <- ext(
      x0,
      min(x0 + tile_ncol * res(r_template)[1], e_aoi$xmax),
      y0,
      min(y0 + tile_nrow * res(r_template)[2], e_aoi$ymax)
    )
  }
}

# --- smoothing + sieve helpers (raster-only) ---
mode3 <- function(x, ...) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_integer_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

smooth3 <- function(x, filename="", overwrite=TRUE) {
  focal(
    x, w=3,
    fun="modal",
    na.rm=TRUE, pad=TRUE,
    filename=filename,
    wopt=list(datatype="INT1U", gdal=c("TILED=YES","COMPRESS=LZW","BIGTIFF=YES")),
    overwrite=overwrite
  )
}

sieve9 <- function(x) {
  p <- patches(x, directions=8)
  f <- as.data.frame(freq(p))
  keep <- f$value[f$count >= 9]
  mask_big <- p %in% keep
  mask(x, mask_big, maskvalues=FALSE)
}

# Run each tile
tile_paths <- character(0)

for (i in seq_along(tile_exts)) {
  ex <- tile_exts[[i]]
  
  # fast reject: if tile doesn't intersect AOI at all
  if (is.null(relate(ext(aoi_v), ex, "intersects"))) {
    # relate() can be picky; a safer check:
    if (!intersect(as.polygons(ex, crs=crs(aoi_v)), aoi_v) |> nrow() > 0) {
      next
    }
  }
  
  message("Tile ", i, "/", length(tile_exts))
  
  preds_tile <- crop(preds, ex)
  preds_tile <- mask(preds_tile, aoi_v)
  
  # if tile is empty after masking, skip
  if (ncell(preds_tile[[1]]) == 0) next
  
  out_raw   <- file.path(out_dir_tiles, sprintf("class_tile_%04d_raw.tif", i))
  out_smooth<- file.path(out_dir_tiles, sprintf("class_tile_%04d_smooth.tif", i))
  out_final <- file.path(out_dir_tiles, sprintf("class_tile_%04d_final.tif", i))
  
  # overwrite to prevent disk bloat
  if (file.exists(out_raw))   file.remove(out_raw)
  if (file.exists(out_smooth))file.remove(out_smooth)
  if (file.exists(out_final)) file.remove(out_final)
  
  # classify tile
  class_raw <- app(
    preds_tile,
    fun=assign_cluster,
    centers=centers,
    mu=mu,
    sd=sdv,
    filename=out_raw,
    wopt=list(datatype="INT1U", gdal=c("TILED=YES","COMPRESS=LZW","BIGTIFF=YES")),
    overwrite=TRUE
  )
  names(class_raw) <- "class_id"
  
  # smooth + sieve (raster)
  class_sm <- smooth3(class_raw, filename=out_smooth, overwrite=TRUE)
  class_fx <- sieve9(class_sm)
  
  writeRaster(
    class_fx, out_final,
    datatype="INT1U",
    gdal=c("TILED=YES","COMPRESS=LZW","BIGTIFF=YES"),
    overwrite=TRUE
  )
  
  # optional: delete intermediates to save space
  if (file.exists(out_raw))   file.remove(out_raw)
  if (file.exists(out_smooth))file.remove(out_smooth)
  
  tile_paths <- c(tile_paths, out_final)
}

stopifnot(length(tile_paths) > 0)

# ============================================================
# 8) Build final mosaic as VRT + optional GeoTIFF
# ============================================================
terra::vrt(tile_paths, filename=out_vrt_class, overwrite=TRUE)

# Optional: materialize to a single GeoTIFF (big!)
# If you only need a mosaic for viewing/processing, keep the VRT.
class_mosaic <- rast(out_vrt_class)
writeRaster(
  class_mosaic, out_tif_class,
  datatype="INT1U",
  gdal=c("TILED=YES","COMPRESS=LZW","BIGTIFF=YES"),
  overwrite=TRUE
)

message("DONE. VRT: ", out_vrt_class)
message("GeoTIFF: ", out_tif_class)
