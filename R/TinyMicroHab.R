# Packages
# Working Directory
library(here)
# Spatial core
library(terra)        # rasters & vectors
library(sf)           # vectors / I/O
# Data
library(dplyr)
# Download & pre-process Sentinel-2 (10 m)
#library(sen2r)        # easiest way to get cloud-masked S2 L2A
# Modeling
library(ranger)       # fast random forest
# Optional: texture features
# install.packages("glcm")
library(glcm)
# Stratified survey design (GRTS)
library(spsurvey)
library(terra)   # v1.7+ recommended
library(stats)
library(here)

#---------------- 0) Inputs ----------------
# naip: SpatRaster (one tile, 4 bands)
naip <- rast(here("data/naip_tiles_2023/az_m_3111038_nw_12_030_20231016_20240119.tif"))

# Find the total rows and columns
nrows_total <- nrow(naip)
ncols_total <- ncol(naip)

# Define the window size
size <- 10000

# Calculate pixel-based coordinates around the center
r_center <- nrows_total / 2
c_center <- ncols_total / 2

# Convert cell indices to spatial coordinates
xy_center <- xyFromCell(naip, cellFromRowCol(naip, r_center, c_center))

# Get the resolution (pixel size)
resx <- res(naip)[1]
resy <- res(naip)[2]

# Build an extent 1000 × 1000 pixels around center
half <- size / 2
ext_small <- ext(
  xy_center[1] - half * resx,
  xy_center[1] + half * resx,
  xy_center[2] - half * resy,
  xy_center[2] + half * resy
)

# Crop that small area
naip_small <- crop(naip, ext_small)

# Plot
plotRGB(naip_small, 1,2,3, stretch="lin")

# aoi : SpatVector or sf object (optional; can skip crop/mask if already clipped)

# DEM 1m
dem <- rast(here("data/dem/USGS_1M_12_x56y349_AZ_CochiseCounty_2020_B20.tif"))

# Crop to same area as NAIP
dem_small <- crop(dem, ext_small)

# Plot
plot(dem_small)
#---------------- 1) Band layout helper ----------------
guess_layout <- function(x) {
  if (nlyr(x) != 4) stop("Expected 4-band NAIP (RGBN). Got ", nlyr(x), " bands.")
  m <- sapply(1:4, function(i) global(x[[i]], "mean", na.rm = TRUE)[[1]])
  nir_idx <- which.max(m)
  if (nir_idx == 4) return("RGBN")
  if (nir_idx == 1) return("NRGB")
  "RGBN"  # fallback
}

layout <- guess_layout(naip_small)  # or set layout <- "RGBN"/"NRGB"
idx <- switch(layout,
              "RGBN" = list(R=1, G=2, B=3, NIR=4),
              "NRGB" = list(NIR=1, R=2, G=3, B=4),
              stop("Unknown layout: ", layout)
)

naip4 <- c(naip_small[[idx$R]], naip_small[[idx$G]], naip_small[[idx$B]], naip_small[[idx$NIR]])
names(naip4) <- c("R","G","B","NIR")

# downsample to 50cm
dem_small <- aggregate(dem_small, fact = 2, fun = mean)

# Ensure CRS for dem and naip match
dem_small <- resample(dem_small, naip4[[1]], method = "bilinear")
dem_small <- mask(crop(dem_small, naip4[[1]]), naip4[[1]])
dem_small <- dem_small[[1]]

#---------------- 2.1) Indices ----------------
ndvi <- (naip4$NIR - naip4$R) / (naip4$NIR + naip4$R)
ndwi <- (naip4$G   - naip4$NIR) / (naip4$G   + naip4$NIR)   # McFeeters (green vs NIR)
savi <- 1.5 * (naip4$NIR - naip4$R) / (naip4$NIR + naip4$R + 0.5)
names(ndvi) <- "NDVI"; names(ndwi) <- "NDWI"; names(savi) <- "SAVI"

# MSAVI2 calcs
# 1) Extract bands (NAIP: 1=Red, 2=Green, 3=Blue, 4=NIR)
red <- naip4[[1]]
nir <- naip4[[4]]

# 2) Convert to reflectance-like scale (MSAVI2 is sensitive to scaling)
redf <- red / 255
nirf <- nir / 255

# 3) MSAVI2
#    MSAVI2 = (2*NIR + 1 - sqrt((2*NIR + 1)^2 - 8*(NIR - RED))) / 2
msavi2 <- (2*nirf + 1 - sqrt((2*nirf + 1)^2 - 8*(nirf - redf))) / 2
msavi2 <- clamp(msavi2, 0, 1)       # keep in [0,1]
names(msavi2) <- "MSAVI2"

# 4) (Optional) local texture of MSAVI2 (5x5 window SD)
#    This helps separate homogeneous vs. heterogeneous vegetation patches
msavi2_sd <- terra::focal(
  msavi2,
  w = matrix(1, 5, 5),
  fun = function(x, ...) sd(x, na.rm = TRUE)
)
names(msavi2_sd) <- "MSAVI2_sd"

# 5) Ensure exact alignment to your template (pick one base, e.g., the DEM or NAIP band 1)
templ <- naip4[[1]]
msavi2    <- resample(msavi2, templ, method = "bilinear")
msavi2_sd <- resample(msavi2_sd, templ, method = "bilinear")


# Texture proxy: focal SD of NDVI (fast). Adjust window if you like.
w <- matrix(1, 5, 5)    # ~1.5 m at 30 cm
ndvi_sd <- terra::focal(
  ndvi,
  w = matrix(1, 5, 5),
  fun = function(x, ...) sd(x, na.rm = TRUE),
  progress = T
)
names(ndvi_sd) <- "NDVI_sd"

#---------------- 2.2) Topographic Derivatives ---------------
slope  <- terrain(dem_small, v = "slope",  unit = "degrees")
aspect <- terrain(dem_small, v = "aspect", unit = "radians")

# Aspect is undefined where slope ~ 0; use a tiny threshold to catch numerics
thresh <- 0.5  # degrees; tweak (0.1–1.0) after you check slope stats

asp_sin <- ifel( is.na(slope) | is.na(aspect) | slope < thresh,
                 0, sin(aspect) )
asp_cos <- ifel( is.na(slope) | is.na(aspect) | slope < thresh,
                 0, cos(aspect) )

# --- 4) Snap to your template AFTERWARDS ----------------------------------
templ   <- naip4[[1]]
slope   <- resample(slope,   templ, method = "bilinear")
asp_sin <- resample(asp_sin, templ, method = "bilinear")
asp_cos <- resample(asp_cos, templ, method = "bilinear")

names(asp_sin) <- "asp_sin"
names(asp_cos) <- "asp_cos"

compareGeom(asp_cos, dem_small, stopOnError = FALSE)

tpi    <- terrain(dem_small, v = "TPI")    # topographic position index
tri    <- terrain(dem_small, v = "TRI")    # terrain ruggedness index

#---------------- 3) Predictor stack & optional clip ----------------
preds <- c(naip4, ndvi, ndwi, ndvi_sd, dem_small, slope, asp_sin, asp_cos, tpi, tri)
# if (exists("aoi")) {
#   aoi_v <- vect(aoi)
#   preds <- mask(crop(preds, aoi_v), aoi_v)
# }

# sanity checks (optional)
# preds; res(preds); nlyr(preds)
# plotRGB(naip4, r=1, g=2, b=3, stretch="lin")
# plot(ndvi)

#---------------- 4) Sample training data ----------------
set.seed(42)
# spatSample returns a SpatRaster with samples; convert to data.frame
samp <- spatSample(preds, size = 50000, method = "random", na.rm = TRUE, as.points = FALSE)
X <- na.omit(as.data.frame(samp, cells = FALSE))  # columns: R,G,B,NIR, NDVI, NDWI, SAVI, NDVI_sd

#---------------- 5) Fit kmeans on scaled features ----------------
k <- 6  # tweak to your expected microhabitats
X_scaled <- scale(X)  # store center/scale for later
sc_center <- attr(X_scaled, "scaled:center")
sc_scale  <- attr(X_scaled, "scaled:scale")

km <- kmeans(X_scaled, centers = k, nstart = 10)
km_centers <- km$centers
mu <- attr(scale(X), "scaled:center")  # use the *same* object you used to fit
sd <- attr(scale(X), "scaled:scale")

#---------------- 6) Classifier for full raster ----------------
# Terra will feed a data.frame (rows=pixels, cols=layers) into this function in chunks.
# We must: (a) apply SAME scaling; (b) assign to nearest km$centers (Euclidean).
# Custom per-pixel assigner: scale inputs -> nearest center
assign_cluster <- function(v, centers, mu, sd) {
  if (any(is.na(v))) return(NA_integer_)
  z <- (v - mu) / sd
  # centers is k x p, z is length p
  d2 <- rowSums((centers - matrix(z, nrow(centers), ncol(centers), byrow=TRUE))^2)
  which.min(d2)
}

# Ensure band order in 'preds' matches the columns of X used to fit kmeans
stopifnot(nlyr(preds) == ncol(km_centers))

# Run prediction chunk-wise; write to disk to avoid RAM spikes
out_tif <- here("output/naip_dem_kmeans6_class.tif")
class_r <- app(
  preds,
  fun = assign_cluster,
  centers = km_centers,
  mu = mu,
  sd = sd,
  filename = out_tif,
  wopt = list(datatype = "INT2U", gdal = c("TILED=YES","COMPRESS=LZW")),
  overwrite = TRUE
)


names(class_r) <- "class_id"

#---------------- 7) Optional: quicklook PNG & a basic palette ----------------
pal <- c("#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a","#b15928")
png(here("figs/naip_dem_kmeans6_quicklook.png"), width = 1400, height = 1200, res = 150)
par(mfrow = c(1, 2))
plot(class_r, col = pal, main = "NAIP/DEM k-means (k=6)")
plotRGB(naip_small, 1,2,3, stretch="lin")
dev.off()

#save a png
# Define your color palette first
pal <- rainbow(6)   # or whatever palette you're using

# Save as PNG
png(here("figs/naip_dem_kmeans6_comparison.png"), width = 2000, height = 1000, res = 150)

par(mfrow = c(1, 2))

# Left: classified raster
plot(class_r, col = pal, main = "NAIP k-means (k=6)")

# Right: true-color NAIP RGB composite
plotRGB(naip_small, r = 1, g = 2, b = 3, stretch = "lin")

dev.off()  # close the PNG device to finalize the file


class_r
class_r <- round(class_r)
terra::writeRaster(class_r, here("output/classified/class_r_INT1U.tif"),
                   datatype = "INT1U",
                   overwrite = TRUE)

#--------------- Smooth and Polygonize -----------------
# Optional: small modal filter to de-speckle (3x3)
terraOptions(todisk = TRUE)  # stream to disk

# 3x3 modal (majority) filter for integer class raster
# - removes NAs from the window
# - returns NA if *all* neighbors are NA
# 0) Define the modal filter safely (returns a single integer)
modal3 <- function(x, ...) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_integer_)
  as.integer(terra::modal(x, ties = "random"))
}

# 1) Make sure class_r is a single-layer integer raster
stopifnot(nlyr(class_r) == 1)
class_r <- as.int(round(class_r))

# 2) (Optional) write your input ONCE (different from output file)
terra::writeRaster(
  class_r,
  here::here("output/classified/class_r_INT1U.tif"),
  datatype = "INT1U",
  overwrite = TRUE
)

# 3) Write focal result to a DIFFERENT filename
out_focal <- here::here("output/classified/class_r_modal3_INT1U.tif")

class_f <- terra::focal(
  class_r,
  w = matrix(1, 3, 3),
  fun = modal3,
  pad = TRUE, padValue = NA,
  filename = out_focal,
  datatype = "INT1U",
  wopt = list(gdal = c("TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES")),
  overwrite = TRUE
)


# Dissolve contiguous pixels into polygons; drop tiny slivers
polys <- as.polygons(class_f, dissolve=TRUE, values=TRUE, na.rm=TRUE)
polys <- st_as_sf(polys) %>% st_make_valid()
polys <- polys %>%
  mutate(area_m2 = as.numeric(st_area(geometry))) %>%
  filter(area_m2 >= 1000)  # keep polygons ≥ 0.1 ha (adjust)
