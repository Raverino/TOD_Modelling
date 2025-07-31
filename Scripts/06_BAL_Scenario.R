# Balanced (BAL) scenario
# This script modifies the Traceback baseline raster to reflect a denser, TOD-aligned urban growth scenario.
# A scenario that is between IB and RETRO

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
TOD_DISTANCE <- 800            # Distance (in meters) to consider as TOD zone
EXPAND_CLASS1 <- 0.15          # % of class 2 pixels converted to class 1
EXPAND_CLASS2 <- 0.25          # % of class 3 pixels converted to class 2
EXPAND_CLASS4 <- 0.15          # % of remaining class 3 to convert to class 4
DECREASE_CLASS5 <- 0.30        # % of Class 5 to convert to Class 7 (green space)
DECREASE_CLASS10 <- 0.15       # % of Class 10 to convert to Class 7 (green space)
CLASS4_BUFFER <- 300           # Maximum distance from class 1 or 2 for class 4 expansion
USE_TOD_DECAY <- TRUE          # If TRUE, apply TOD influence as decaying function
TOD_DECAY_RATE <- 0.005        # Decay rate of TOD influence with distance
SCENARIO_NAME <- "BAL"         # Scenario tag for output naming
DOWNSAMPLE_FACTOR <- 4         # Factor to downsample for distance calculation

# Extract input name from loaded raster
input_name <- tools::file_path_sans_ext(basename(urban_path))

# Copy baseline
urban_BAL <- urban_baseline

# Define masks
class1_mask <- urban_BAL == 1
class2_mask <- urban_BAL == 2
class3_mask <- urban_BAL == 3
class4_mask <- urban_BAL == 4
class5_mask <- urban_BAL == 5
class6_mask <- urban_BAL == 6
class7_mask <- urban_BAL == 7
class8_mask <- urban_BAL == 8
class9_mask <- urban_BAL == 9
class10_mask <- urban_BAL == 10

# Define TOD zone 
tod_zone <- dist_transit <= TOD_DISTANCE

#  Zoning Preparation 
model_zoning <- zoning
model_zoning[is.na(model_zoning)] <- 1  # Assume NA = developable

# Ensure existing urban (classes 1–6) - Classes that are going to expand
existing_urban <- urban_baseline
existing_urban[!(existing_urban[] %in% 1:6)] <- 0
existing_urban[existing_urban > 0] <- 1
model_zoning <- max(model_zoning, existing_urban, na.rm = TRUE)

# Transformations 
set.seed(42)

# 1. Class 1 grows by % into class 2 areas (TOD + zoning)
tod_class2 <- which((urban_BAL[] == 2) & tod_zone[] & model_zoning[] == 1)
n_class1_expand <- floor(length(tod_class2) * EXPAND_CLASS1)
expanded_class1 <- sample(tod_class2, n_class1_expand)
ur_class2_for_1 <- expanded_class1
urban_BAL[ur_class2_for_1] <- 1

# 2. Class 2 grows by % into class 3 areas (TOD + zoning)
tod_class3 <- which((urban_BAL[] == 3) & tod_zone[] & model_zoning[] == 1)
n_class2_expand <- floor(length(tod_class3) * EXPAND_CLASS2)
expanded_class2 <- sample(tod_class3, n_class2_expand)
ur_class3_for_2 <- expanded_class2
urban_BAL[ur_class3_for_2] <- 2

# 3. Class 4 takes % of remaining Class 3 (zoning + within buffer of class 1/2) using optimized approach
urban_lowres <- aggregate(urban_BAL, fact=DOWNSAMPLE_FACTOR, fun=modal)
ur_class_mask_lowres <- (urban_lowres == 1 | urban_lowres == 2)
ur_class_dist_lowres <- distance(ur_class_mask_lowres)
ur_class_dist <- resample(ur_class_dist_lowres, urban_BAL, method="bilinear")
buffer_zone <- ur_class_dist <= CLASS4_BUFFER
remaining_class3 <- which((urban_BAL[] == 3) & model_zoning[] == 1 & buffer_zone[] == 1)
n_class4_replace <- floor(length(remaining_class3) * EXPAND_CLASS4)
ur_class3_for_4 <- sample(remaining_class3, n_class4_replace)
urban_BAL[ur_class3_for_4] <- 4

# 4. Class 5 decrease by 15% to Class 7 (green space conversion)
class5_in_tod <- which((urban_BAL[] == 5) & tod_zone[] & model_zoning[] == 1)
n_class5_decrease <- floor(length(class5_in_tod) * DECREASE_CLASS5)
selected_class5 <- sample(class5_in_tod, n_class5_decrease)
urban_BAL[selected_class5] <- 7

# 5. Class 10 decrease by 40% to Class 7 (green space conversion)
class10_in_tod <- which((urban_BAL[] == 10) & tod_zone[] & model_zoning[] == 1)
n_class10_decrease <- floor(length(class10_in_tod) * DECREASE_CLASS10)
selected_class10 <- sample(class10_in_tod, n_class10_decrease)
urban_BAL[selected_class10] <- 7


# Save Result 
output_dir <- "C:/Users/darth/OneDrive/Dokumente/Pessoal/UTMaster/MasterThesis/QGIS_MasterThesis/tif/Aligned_Apolus/BAL_Scenario/"
writeRaster(urban_BAL, filename=paste0(output_dir, input_name, "_", SCENARIO_NAME, ".tif"), format="GTiff", overwrite=TRUE)