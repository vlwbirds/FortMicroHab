# Packages
library(terra)   # v1.7+ recommended
library(stats)
library(mapview)
library(here)

#---------------- 0) Inputs ----------------
# naip: SpatRaster (one tile, 4 bands)
naip <- rast("data/naip_tiles_2023/az_m_3111038_nw_12_030_20231016_20240119.tif")

# list subdatasets if there are multiple
info <- terra::sds(naip)
info

# img_1 <- rast(info[1])
# img_2 <- rast(info[2])

# Find the total rows and columns
nrows_total <- nrow(naip)
ncols_total <- ncol(naip)

# Define the window size
size <- 1000

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
dem <- rast("data/dem/USGS_1M_12_x56y349_AZ_CochiseCounty_2020_B20.tif")

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

# calculate naip reflectances
reflect_naip <- naip4 / 255

# downsample to 50cm
dem_small <- aggregate(dem_small, fact = 2, fun = mean)

# Ensure CRS for dem and naip match
dem_small <- resample(dem_small, naip4[[1]], method = "bilinear")
dem_small <- mask(crop(dem_small, naip4[[1]]), naip4[[1]])
dem_small <- dem_small[[1]]

#------------- 2) Calculate SVI (shaded vegetation index) & RVI (Ratio Vegetation Index) ----
rvi <- app(reflect_naip, fun = function(x) {
  nir <- x[4]   # assuming band order is R,G,B,NIR — adjust if different
  red <- x[1]
  ifelse(red == 0 | is.na(red), NA, nir / red)
})
names(rvi) <- "RVI"

svi <- app(reflect_naip, fun = function(x) {
  red <- x[1]
  ifelse(red == 0 | is.na(red), NA, 1 / red)
})
names(svi) <- "SVI"

par(mfrow = c(1,2))
hist(svi)
hist(rvi)
