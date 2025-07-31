# This script allows to extract basic statistics that relate each urban pixel to the urban area situated
# It serves as a base for the Traceback Algorithm

library(raster)
library(ggplot2)

setwd("C:/Users/darth/OneDrive/Dokumente/Pessoal/UTMaster/MasterThesis/QGIS_MasterThesis/tif/aligned_apolus")

# Load your most recent land use map
landuse_map <- raster("Roma_capitale_2018.tif")

# Load spatial criteria maps
transit_dist <- raster("TransitNetwork_Aligned.tif")   # Distance to transit network (euclidean)
road_dist <- raster("RoadNetwork_Aligned.tif")         # Distance to road network (euclidean)
slope <- raster("Slope_Aligned.tif")                   # Slope map (degrees)
city_center_dist <- raster("urban_polos_distance.tif") # Distance to the urban center (euclidean)

# Check if files exist
if (!file.exists("Roma_2018_11class.tif")) stop("❌ Land use raster file not found.")
if (!file.exists("TransitNetwork_Aligned.tif")) stop("❌ Transit raster file not found.")
if (!file.exists("RoadNetwork_Aligned.tif")) stop("❌ Road raster file not found.")
if (!file.exists("Slope_Aligned.tif")) stop("❌ Slope raster file not found.")
if (!file.exists("urban_polos_distance.tif")) stop("❌ Urban centers distance file not found.")

# Function to analyze spatial characteristics by land use class
analyze_landuse_characteristics <- function(landuse_map, criteria_list) {
  # Get unique land use classes
  classes <- unique(getValues(landuse_map))
  classes <- classes[!is.na(classes)]
  
  # Create empty results list
  results <- list()
  
  # For each land use class
  for (class_id in classes) {
    cat("Processing class", class_id, "\n")
    # Create mask for this class
    class_mask <- landuse_map == class_id
    
    # Analyze each criteria for this class
    class_stats <- list()
    for (criteria_name in names(criteria_list)) {
      criteria_raster <- criteria_list[[criteria_name]]
      
      # Extract values where this class exists
      # Get cells where mask is TRUE (value is 1)
      cells <- which(values(class_mask) == 1)
      
      # If we have cells with this land use class
      if (length(cells) > 0) {
        # Extract criteria values at these cell positions
        values <- values(criteria_raster)[cells]
        values <- values[!is.na(values)]
        
        # Calculate statistics if we have values
        if (length(values) > 0) {
          class_stats[[criteria_name]] <- list(
            mean = mean(values, na.rm=TRUE),
            median = median(values, na.rm=TRUE),
            q25 = quantile(values, 0.25, na.rm=TRUE),
            q75 = quantile(values, 0.75, na.rm=TRUE),
            min = min(values, na.rm=TRUE),
            max = max(values, na.rm=TRUE)
          )
        } else {
          class_stats[[criteria_name]] <- list(
            mean = NA, median = NA, q25 = NA, q75 = NA, min = NA, max = NA
          )
          warning(paste("No valid values found for class", class_id, "and criteria", criteria_name))
        }
      } else {
        class_stats[[criteria_name]] <- list(
          mean = NA, median = NA, q25 = NA, q75 = NA, min = NA, max = NA
        )
        warning(paste("No cells found for class", class_id))
      }
    }
    
    results[[as.character(class_id)]] <- class_stats
  }
  
  return(results)
}

# Create list of criteria rasters
criteria <- list(
  transit_dist = transit_dist,
  road_dist = road_dist,
  slope = slope,
  city_center_dist = city_center_dist
)

# Print some basic information about our rasters
print("Checking raster dimensions and CRS compatibility:")
print(paste("Land use map dimensions:", nrow(landuse_map), "x", ncol(landuse_map)))
for (name in names(criteria)) {
  print(paste(name, "dimensions:", nrow(criteria[[name]]), "x", ncol(criteria[[name]])))
  print(paste(name, "CRS matches land use map:", identical(crs(criteria[[name]]), crs(landuse_map))))
}

# Run analysis
landuse_characteristics <- analyze_landuse_characteristics(landuse_map, criteria)

# Print summary statistics for each class
for (class_id in names(landuse_characteristics)) {
  cat("====== Land Use Class:", class_id, "======\n")
  for (criteria_name in names(landuse_characteristics[[class_id]])) {
    cat(criteria_name, ":\n")
    stats <- landuse_characteristics[[class_id]][[criteria_name]]
    cat("  Mean:", stats$mean, "\n")
    cat("  Median:", stats$median, "\n")
    cat("  25% Quantile:", stats$q25, "\n")
    cat("  75% Quantile:", stats$q75, "\n")
    cat("  Range:", stats$min, "-", stats$max, "\n\n")
  }
}

# Convert nested list to data frame
stats_df <- do.call(rbind, lapply(names(landuse_characteristics), function(class_id) {
  class_stats <- landuse_characteristics[[class_id]]
  do.call(rbind, lapply(names(class_stats), function(criteria_name) {
    stats <- class_stats[[criteria_name]]
    data.frame(
      class_id = class_id,
      criteria = criteria_name,
      mean = stats$mean,
      median = stats$median,
      q25 = stats$q25,
      q75 = stats$q75,
      min = stats$min,
      max = stats$max
    )
  }))
}))

# Save to CSV file
write.csv(stats_df, "landuse_statistics.csv", row.names = FALSE)
cat("✅ The file has been created.\n")