#Traceback Algorithm
#This script allows to transform a binary class land-use map into multiple classes based on previous land-use statistics

library(raster)

cat("📁 Setting working directory for rasters...\n")
setwd("C:/Users/darth/OneDrive/Dokumente/Pessoal/UTMaster/MasterThesis/QGIS_MasterThesis/tif/halfway")

# Load simulated urban/non-urban map
urban_file <- "TM_S_R_No_Zoning.tif"
if (!file.exists(urban_file)) stop("❌ Urban raster file not found.")
urban_2050 <- raster(urban_file)
cat("✅ Urban raster loaded successfully.\n")

# Load spatial criteria maps
transit_file <- "Transit_Network.tif"
road_file <- "Road_Network.tif"
slope_file <- "Slope.tif"
city_file <- "Urban_Polos.tif"

if (!file.exists(transit_file)) stop("❌ Transit raster file not found.")
if (!file.exists(road_file)) stop("❌ Road raster file not found.")
if (!file.exists(slope_file)) stop("❌ Slope raster file not found.")
if (!file.exists(city_file)) stop("❌ Urban centers distance file not found.")

transit_dist <- raster(transit_file)
road_dist <- raster(road_file)
slope <- raster(slope_file)
urban_dist <- raster(city_file)

cat("✅ All spatial criteria rasters loaded successfully.\n")

# Check raster alignment
if (!compareRaster(urban_2050, transit_dist, road_dist, slope, urban_dist, stopiffalse = FALSE)) {
  stop("❌ Rasters are not aligned or have mismatched resolution/extents.")
}
cat("✅ Raster alignment check passed.\n")

# Extract values from rasters
urban_vals   <- getValues(urban_2050)
transit_vals <- getValues(transit_dist)
road_vals    <- getValues(road_dist)
slope_vals   <- getValues(slope)
urban_dist <- getValues(urban_dist)

# Initialize classified values
classified_vals <- rep(NA, length(urban_vals))

# Define urban cells
urban_cells <- which(urban_vals == 1)
cat(paste0("🔍 Found ", length(urban_cells), " urban cells to classify.\n"))

# Spatial criteria for urban cells
dist_t <- transit_vals[urban_cells]  #Distance to transit
dist_r <- road_vals[urban_cells]     #Distance to road network
dist_u <- urban_dist[urban_cells]    #Distance to urban polos
slope_u <- slope_vals[urban_cells]   #Distance to slope


# The magic starts here...
class <- rep(NA, length(urban_cells))
# Classification rules based on statistical analysis
class[is.na(class) & dist_t <= 1160 & dist_r <= 30 & slope_u <= 6.5 & dist_u <= 80] <- 1      # High-Density
class[is.na(class) & dist_t <= 3000 & dist_r <= 42 & slope_u <= 5.7 & dist_u <= 120] <- 2     # Medium-Density
class[is.na(class) & dist_t <= 3880 & dist_r <= 60 & slope_u <= 7.2 & dist_u <= 240] <- 3     # Low-Density
#class[is.na(class) & dist_t > 3700 & dist_r <= 86 & dist_u > 200] <- 4                        # Old Isolated Structures
class[is.na(class) & dist_r <= 85 & dist_t <= 2700 & slope_u <= 5.0 & dist_u <= 222] <- 4     # Commercial/Industrial
class[is.na(class) & dist_r <= 10] <- 5                                                       # Transportation
class[is.na(class) & dist_t <= 2900 & dist_r <= 72 & slope_u <= 5.1 & dist_u <= 172] <- 6     # Dev Zones
class[is.na(class) & dist_t <= 2000 & dist_r > 50 & dist_u <= 150] <- 7                       # Urban Green
class[is.na(class) & dist_t > 3000 & dist_r > 100 & dist_u > 175] <- 10                       # Excluded Zones
# Handle remaining NA using proximity to existing classes
# Handle any unclassified urban cells
class[is.na(class) & dist_r <=100] <- 4
class[is.na(class)] <- 2


# And then for the non-urban cells we got the classification
# Define non-urban cells
non_urban_cells <- which(urban_vals == 0)

# Extract spatial values for non-urban cells
non_urban_dist_t <- transit_vals[non_urban_cells]
non_urban_dist_r <- road_vals[non_urban_cells]
non_urban_slope  <- slope_vals[non_urban_cells]
non_urban_dist_u <- urban_dist[non_urban_cells]

# Initialize classification vector
non_urban_class <- rep(NA, length(non_urban_cells))

non_urban_class[non_urban_slope <= 6.5 & non_urban_dist_r <= 212] <- 8 # Class 8: Agricultural & Rural Areas
non_urban_class[is.na(non_urban_class)] <- 9 # Class 9: Water & Natural Areas (fallback)

# Assign to classified raster
classified_vals[non_urban_cells] <- non_urban_class


# Assign classes to raster
classified_vals[urban_cells] <- class
detailed_classes <- urban_2050
values(detailed_classes) <- classified_vals

cat("✅ Classification complete.\n")

# Save output
output_file <- sub("\\.tif$", "_TB.tif", urban_file)
writeRaster(detailed_classes, output_file, overwrite = TRUE)

# Confirm file creation
if (file.exists(output_file)) {
  cat(paste0("🎉 Raster successfully created: ", output_file, "\n"))
} else {
  cat("❌ Raster write failed.\n")
}
