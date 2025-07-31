library(raster)
setwd("C:/Users/darth/OneDrive/Dokumente/Pessoal/UTMaster/MasterThesis/QGIS_MasterThesis/tif/aligned_apolus/BAL_scenario")

# Load map for urban pixels output
map_selected <- brick("TB_R_Zoning_BAL.tif")

# Extract the land-use band (assumed to be the first band)
landuse <- map_selected[[1]]

# Function to count the pixels for each land-use class
count_pixels <- function(raster_layer) {
  pixel_counts <- table(values(raster_layer))
  return(pixel_counts)
}

# Count pixels for the map selected
pixel_counts_landuse <- count_pixels(landuse)

# Print the pixel counts for each year
cat("Pixel counts in the map:\n")
print(pixel_counts_landuse)


#Create a csv file with the data in horizontal
write.csv(t(as.data.frame(pixel_counts_landuse)), file = "pixel_counts_landuse.csv", row.names = FALSE)
landuse_csv <- ("pixel_counts_landuse.csv")

# Confirm file creation
if (file.exists(landuse_csv)) {
  cat(paste0("🎉 The data is now saved and available: ", landuse_csv, "\n"))
} else {
  cat("❌ RTRy again, something wrong happened\n")
}