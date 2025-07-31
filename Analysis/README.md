
# Analysis

This folder contains processed data and summary tables used to evaluate the outcomes of the TOD simulation framework. The primary focus is on quantifying housing potential and carbon emissions under each scenario. Three files can be found in these folders, **one being the core** for the analysis and **two are support files**. The description on how these work can be found below.

---

## 📁 Contents

### `Data_Analysis.xlsx` ✅ **Core file**
This is the main results table and forms the analytical core of the thesis.
The input data comes from GIS analysis and literature review. The cells in white are filled automatically, so do only change those in orange (from GIS) or green (literature review). 
For the **urban pixels** (U.P) it has been used the script PixelCounter for each land-use. For **TOD U.P.** and **Housing Capacity** the method of Zonal Histogram in GIS has been used. 

The method for **TOD U.P.** was the following:
1) Creation of a vector layer that includes the public transport (BUS, tram and train stations... depends on the detail of the analysis).
2) A buffer of 1000m (if decided to go in more detail, read TOD requirements). This represents the area of TOD. 
3) Zonal histogram between residential areas (class 1,2 and 3) layer and the TOD buffer. This allows to tell how much people are served by TOD in this area.

As for the **Housing capacity** was the following:
1) OSM query for building (apartments and house type)
2) Zonal histogram relation between these building types and the different 3 classes. 

**Sheets include:**
- **Transit Accessibility**: It includes Transit Urban Coverage in general & within TOD distance (1000m); Population Served; The possible car trips reduced and their respective carbon emissions reduced.
- **Housing Capacity** per scenario (based on reclassified land-use pixels × density factors)
- **Urban Carbon Emissions** estimates based on land-use type and extent
- **Carbon Emissions Reduction Potential**, comparing baseline and TOD scenarios

This file directly supports the findings in **Chapter 5 (Results)** and **Chapter 6 (Discussion)** of the thesis.

---

### `Urban_Pixel_Data.xlsx` 📐 *(Auxiliary)*
Tracks how many pixels each map contains per land-use class and scenario, and converts that into:
- Total area (in km²)
- Percentage change between scenarios
- Zoning Percentage

Used to support spatial consistency checks and validate simulation coverage.

---

### `Landuse_Stats.xlsx` 🔁 *(Auxiliary)*
An example of the Stats_Raster output used for the **traceback algorithm**.

Includes:
- Tabulated suitability statistics (distance to features, slope)
---

## Citation

If you use this data in your own academic work, please cite:

Paulo Dimas. (2025). TOD Simulation Scripts and Data (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.16635783
