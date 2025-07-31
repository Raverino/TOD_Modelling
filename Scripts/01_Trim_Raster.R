#This script allows to trim a raster to it's extent

library(terra)
# Load your raster
setwd("C:/Users/darth/OneDrive/Dokumente/Pessoal/UTMaster/MasterThesis/QGIS_MasterThesis/tif/aligned_apolus")
input_filename <- "TM_S_60.tif"
r <- rast(input_filename)

# Create a mask where values are not black (not equal to 0)
mask <- r > 0

# Get the extent of non-zero values
ext <- ext(trim(mask))

# Crop the original raster to this extent
result <- crop(r, ext)

# Write the result to a new file
writeRaster(result, sub("\\.tif$", "_trimmed.tif", input_filename), overwrite=TRUE)
