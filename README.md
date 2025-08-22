# TOD-Simulation: Transit-Oriented Development Land-Use Scenario Modeling

This repository contains spatial simulation scripts and supporting files developed for the master's thesis:

**"Modeling Transit-Oriented Development: A Land-Use Approach for Urban Decarbonization"**  
*Paulo Dimas, 2025*
Thesis available at: https://purl.utwente.nl/essays/108255
Master’s Thesis, University of Twente

The simulation framework adapts and extends the SIMLANDER land-use change model using a set of original R scripts to evaluate decarbonization scenarios based on Transit-Oriented Development (TOD) principles. The main case study is the metropolitan area of Rome.

---

## 🗂️ Repository Structure

```
/
├── LICENSE              ← MIT license for all original code
├── README.md            ← Project overview (this file)
├── /scripts/            ← All R scripts used for modeling and scenario simulation
├── /analysis/           ← All Excel Sheets used for the analysis
├── /simlander_mod/      ← Modified script from SIMLANDER.
```
## Online Files

Due to the large size of the examples, you can find it online through the [link](https://drive.google.com/drive/folders/1bI6lVB7FmtDeNoXZKM9h4__mRc65CiLl?usp=sharing) 
```
├── /inputs/             ← Land-use maps, accessibility layers, zoning constraints
├── /outputs/            ← Simulation results, emissions tables, final rasters
---
```
## 📜 Description of Scripts

This repository includes eight R scripts (three of them are scenarios) organized to perform land-use simulation modeling, accessibility evaluation, and TOD scenario execution.

### 1. `Trim_Raster.R`
Trims all input rasters to a consistent spatial extent and resolution.  
Ensures slope, road, transit, and land-use rasters are spatially aligned.

### 2. `Stats_Raster.R`
Generates a `.csv` file with per-pixel spatial statistics. As input it should be the latest land-use map. The features are:
- Transit network
- Road network
- Urban center
- Slope

### 3. `TB_Vector.R` (Traceback Algorithm)
Converts binary urban/non-urban maps into multi-class land-use rasters.  Uses output from `Stats_Raster.R`. 

### 4. Scenario Scripts
- `Scenario_IB.R`: Applies TOD densification logic (Intense Building)
- `Scenario_RETRO.R`: Implements green retrofitting strategy
- `Scenario_BAL.R`: Balanced strategy combining TOD and retrofitting

Each scenario applies rule-based transformations based on spatial suitability metrics and zoning constraints.

### 5. `PixelCounter.R`
Counts the number of pixels per land-use class.  
Useful both before and after scenario simulation to assess land-use shifts.

---

## 🔁 Modeling Workflow

1. Prepare raster data (`Trim_Raster.R`)
	1. Land-use maps with 2 classes (1 = urban and 0 = non-urban)
	2. Road network (euclidean distance)
	3. Transit network (euclidean distance)
	4. Urban distance (euclidean distance)
	5. Slope terrain (degrees)
	6. Zoning/Protected areas with 2 classes (0 = protected and 1 = developable)
2. Compute spatial metrics (`Stats_Raster.R`)
3. Apply traceback to reclassify urban areas (`TB_Vector.R`)
4. Run scenario of choice (IB, RETRO, BAL)
5. Count results (`PixelCounter.R`)
6. Analyse through the different Excel Sheets and GIS software. 

---

## 📍 Data Requirements

All input rasters must be in **the same projection, resolution and aligned**. 

---

## 📄 License

This repository is licensed under the [MIT License](LICENSE).

> This work includes adapted logic inspired by the SIMLANDER land-use simulation model developed by Richard Hewitt, which is described by the author as Free and Open Source Software (FOSS).  
> The SIMLANDER source code itself is not included in this repository. Only original or adapted components are shared here under MIT terms.

---

## 📚 Citation

If you use this code in academic research, please cite:

Paulo Dimas. (2025). TOD Simulation Scripts and Data (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.16635783

**BibTeX**
```bibtex
@mastersthesis{dimas2025tod,
  title     = {Modeling Transit-Oriented Development: A Land-Use Approach for Urban Decarbonization},
  author    = {Paulo Dimas},
  year      = {2025},
  school    = {University of Twente},
  url       = {https://github.com/Raverino/TOD_modelling}
}
```

And optionally cite the SIMLANDER model:

> Hewitt, R. (2022). SIMLANDER: A land-use simulation model. *International Journal of Geographical Information Science*. https://doi.org/10.1080/13658816.2022.2098299

---

## 📬 Contact

📧 For questions, feedback, or collaboration: p.r.almeidadimas@student.utwente.nl


