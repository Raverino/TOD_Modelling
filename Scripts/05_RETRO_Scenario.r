# Retrofitting Scenario (RETRO)
# This script modifies the Traceback baseline raster to reflect the RETRO scenario, 
# Focusing on reuse and transformation of underutilized land within TOD zones.

library(raster)
library(sp)
library(RColorBrewer)

# Setup
# Set working directory to where your raster files are stored
setwd("C:/Users/darth/OneDrive/Dokumente/Pessoal/UTMaster/MasterThesis/QGIS_MasterThesis/tif/halfway")

# Load Raster Data
urban_path <- "TB_R_Zoning.tif"
urban_baseline <- raster(urban_path)
dist_transit <- raster("Transit_Network.tif")  # in meters
zoning <- raster("zoning_rome.tif")     # 0 = protected, 1 = developable

# Parameters
TOD_DISTANCE <- 1000           # Distance (in meters) to consider as TOD zone
REPURPOSE_CLASS4 <- 0.30       # % of Class 4 to convert to Class 1 or 2
REPURPOSE_CLASS10 <- 0.50      # % of Class 10 to convert to Class 7 (green) or 1
REPURPOSE_CLASS5 <- 0.20       # % of Class 5 (infrastructure) to Class 7 (green)
USE_TOD_DECAY <- TRUE          # If TRUE, apply TOD influence as decaying function
TOD_DECAY_RATE <- 0.005        # Decay rate of TOD influence with distance
SCENARIO_NAME <- "RETRO"       # Scenario tag for output naming

# Extract input name from loaded raster
input_name <- tools::file_path_sans_ext(basename(urban_path))

# Copy baseline
urban_RETRO <- urban_baseline

# === Define masks ===
class1_mask <- urban_RETRO == 1
class2_mask <- urban_RETRO == 2
class3_mask <- urban_RETRO == 3
class4_mask <- urban_RETRO == 4
class5_mask <- urban_RETRO == 5
class6_mask <- urban_RETRO == 6
class7_mask <- urban_RETRO == 7
class8_mask <- urban_RETRO == 8
class9_mask <- urban_RETRO == 9
class10_mask <- urban_RETRO == 10

# Define TOD zone
tod_zone <- dist_transit <= TOD_DISTANCE

# Zoning Preparation
model_zoning <- zoning
model_zoning[is.na(model_zoning)] <- 1  # Assume NA = developable

# Ensure existing urban (classes 1–6) - Classes that are going to expand
existing_urban <- urban_baseline
existing_urban[!(existing_urban[] %in% 1:6)] <- 0
existing_urban[existing_urban > 0] <- 1
model_zoning <- max(model_zoning, existing_urban, na.rm = TRUE)

# Transformations
set.seed(42)

# 1. Class 4 (Commercial & Industrial) to Class 1 and 2 (TOD + zoning)
tod_class4 <- which((urban_RETRO[] == 4) & tod_zone[] & model_zoning[] == 1)
n_class4_change <- floor(length(tod_class4) * REPURPOSE_CLASS4)
if(n_class4_change > 0) {
  selected_class4 <- sample(tod_class4, n_class4_change)
  half_point <- floor(n_class4_change/2)
  # Half to residential (class 1)
  urban_RETRO[selected_class4[1:half_point]] <- 1
  # Half to mixed use (class 2)
  urban_RETRO[selected_class4[(half_point+1):n_class4_change]] <- 2
}

# 2. Class 10 (Excluded Industrial & Special) to Class 1 and 7 (TOD + zoning)
tod_class10 <- which((urban_RETRO[] == 10) & tod_zone[] & model_zoning[] == 1)
n_class10_change <- floor(length(tod_class10) * REPURPOSE_CLASS10)
if(n_class10_change > 0) {
  selected_class10 <- sample(tod_class10, n_class10_change)
  half_point <- floor(n_class10_change/2)
  # Half to green space (class 7)
  urban_RETRO[selected_class10[1:half_point]] <- 7
  # Half to residential (class 1)
  urban_RETRO[selected_class10[(half_point+1):n_class10_change]] <- 1
}

# 3. Class 5 (Infrastructure) to Class 7 (Green Space) (TOD + zoning)
tod_class5 <- which((urban_RETRO[] == 5) & tod_zone[] & model_zoning[] == 1)
n_class5_change <- floor(length(tod_class5) * REPURPOSE_CLASS5)
if(n_class5_change > 0) {
  selected_class5 <- sample(tod_class5, n_class5_change)
  urban_RETRO[selected_class5] <- 7
}

# Save Result
output_dir <- "C:/Users/darth/OneDrive/Dokumente/Pessoal/UTMaster/MasterThesis/QGIS_MasterThesis/tif/Aligned_Apolus/Retro_Scenario/"
writeRaster(urban_RETRO, filename=paste0(output_dir, input_name, "_", SCENARIO_NAME, ".tif"), format="GTiff", overwrite=TRUE)