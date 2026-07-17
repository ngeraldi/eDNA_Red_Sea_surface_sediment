# eDNA Red Sea Surface Sediment

R-based eDNA metabarcoding analysis of surface sediments from Red Sea coastal habitats (coral, seagrass, mangrove), sampled in 2016–2017. Two 18S primers are used: **18smini** (18S V9) and **euka02** (18S V7).

All active scripts are in `2026 RS eDNA/`. Run them in order as described below.

---

## Scripts

### Step 1 — `Red_Sea_post_DADA2_filter_and_tidy_withdecon.Rmd`
Imports raw DADA2 outputs, clusters ASVs at 97%, filters by occurrence, blank contamination, and read proportions, rarefies, and exports the main filtered dataset (`RS_surface_filtered_data_all.csv`).

### Step 2 — `Red_sea_external_data.Rmd`
Extracts environmental covariates (SST, nitrate, salinity, primary productivity, distance to shore) from rasters using `terra`/`sf`, and fetches species traits from WoRMS via `worrms`. Exports `Red_sea_eDNA_enviro_var.csv`.

### Step 3 — `Red_sea_analysis_plots_euk.Rmd`
Main community analysis script. Runs NMDS ordination, linear mixed models, spatial autocorrelation, multivariate GLM (`manyglm`), RLQ/fourth-corner trait analysis, and distance decay. Produces all manuscript figures and tables.

### `worms_name_check.R`
Standalone script that validates species names against WoRMS accepted names in batches. Requires internet. Run from the R console after step 1.

---

## Data files

| File | Description |
|------|-------------|
| `RS_surface_sample_data_all.xlsx` | Sample metadata (coordinates, habitat, season, depth); one sheet per primer plus a `sample_data` sheet |
| `RS_surface_filtered_data_all.csv` | Main filtered and rarefied dataset produced by step 1; primary input for step 3 |
| `Red_sea_eDNA_enviro_var.csv` | Environmental covariates per sample produced by step 2 |
| `Red_sea_eDNA_traits_cleaned.csv` | Ecological trait data (trophic level, zone, mobility, larval form) for detected taxa; ecological corrections applied and documented separately |
| `through_water_dist_km.csv` | Pairwise through-water distances (km) between sample sites, computed via NOAA ETOPO1 bathymetry with land as a barrier; used in the distance decay analysis |

---

## Requirements

- R ≥ 4.3
- Key packages: `tidyverse`, `vegan`, `ade4`, `mvabund`, `lme4`, `lmerTest`, `spdep`, `spatialreg`, `mgcv`, `terra`, `sf`, `marmap`, `worrms`, `ggplot2`, `patchwork`, `wesanderson`, `flextable`

---

## Reference

Geraldi et al. (in prep.) — Red Sea coastal sediment eDNA metabarcoding.

For questions contact: nategeraldi@gmail.com
