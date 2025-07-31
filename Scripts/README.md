# Scripts

This folder contains all R scripts used to process raster inputs, generate spatial suitability statistics, reclassify land-use classes, and run TOD-based scenario simulations. The codebase builds on the SIMLANDER model structure and introduces several adaptations for Transit-Oriented Development (TOD) modeling.

All scripts are modular and documented to allow sequential execution or independent testing.

---

## 🗂️ Script Overview

### `01_Trim_Raster.R`
Trims all raster layers (land use, slope, road, transit, zoning) to a shared extent and resolution.  
**Output:** Consistently aligned input rasters used for all further steps.

---

### `02_Stats_Raster.R`
Calculates pixel-level spatial statistics to feed into the scenario algorithms.  
Metrics include:
- Distance to transit
- Distance to roads
- Distance to center
- Slope category

**Output:** `.csv` with per-pixel feature scores used to guide land-use change.

---

### `03_TB_Vector.R` (Traceback Algorithm)
Reclassifies urban pixels in a binary land-use raster into more detailed land-use categories.  
Used after simulation to redistribute land-use types based on spatial suitability.

**Input:** Binary raster (urban / non-urban)  
**Output:** Multi-class raster with classes from `landuse_class_legend.txt`

---

### `04_Scenario_IB.R`
TOD Intense Building scenario.  
Prioritizes high-density development around transit hubs.

### `05_Scenario_RETRO.R`
Green Retrofitting scenario.  
Converts roads and low-density areas into green infrastructure and housing reuse.

### `06_Scenario_BAL.R`
Balanced scenario combining TOD intensification with retrofitting strategy.

All three scenario scripts apply:
- Suitability-based selection using outputs from `Stats_Raster.R`
- Constraints from zoning masks
- Reallocation rules per scenario

**Output:** Future land-use raster after scenario transformation

---

### `07_PixelCounter.R`
Counts the number of pixels per land-use class in any raster file.  
Used before/after scenario simulation to evaluate growth, densification, and conversion.

**Output:** Tabulated summary (`.csv`) of land-use class frequencies.

---

## ⚠️ Dependencies

- **R packages required**:
  - `terra`
  - `raster`
  - `sf`
  - `tidyverse`
  - `data.table`

Ensure all input rasters are in **EPSG:3035** and **100m resolution**.

---

## 🔁 Suggested Execution Order

1. `01_Trim_Raster.R`
2. `02_Stats_Raster.R`
3. `03_TB_Vector.R`
4. One of:
   - `04_Scenario_IB.R`
   - `05_Scenario_RETRO.R`
   - `06_Scenario_BAL.R`
5. `07_PixelCounter.R`

Each script is documented with comments to allow reuse or adaptation to other cities.

---

## 📚 Citation

If you use these scripts or any adaptation of them, please cite:

> Dimas, P. (2025). *Modeling Transit-Oriented Development: A Land-Use Approach for Urban Decarbonization*. University of Twente. 
