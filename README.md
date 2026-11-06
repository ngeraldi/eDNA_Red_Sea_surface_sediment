# eDNA Red Sea Surface Sediment

R-based eDNA metabarcoding analysis of surface sediments from Red Sea coastal habitats (coral, seagrass, mangrove), sampled in 2016–2017. Two 18S primers are used: **18smini** (18S V9) and **euka02** (18S V7).

All active scripts are in `2026 RS eDNA/`. Run them in order as described below.

---

## Scripts

### Step 1 — `Red_Sea_post_DADA2_filter_and_tidy_withdecon.Rmd`
Imports DADA2 outputs and clusters ASVs at 97% (`DECIPHER`). Filters 97% clusters by occurrence (≥ 3 samples) and read share (≥ 0.005% of sample reads), keeps ASVs named at class level or below, removes contaminants (`decontam` plus a blank rule) and non-marine taxa (rules in `Red_sea_non_marine_taxa.csv`), and drops samples with < 1,000 reads in either primer. Exports the filtered dataset (`RS_surface_filtered_data_all.csv.gz`), the dropped-sample table and HTML filter summaries (not committed; regenerated on knitting). Also holds the extraction-replicate analysis.

### Step 2 — `Red_sea_external_data.Rmd`
Extracts environmental covariates (SST, nitrate, salinity, primary productivity, distance to shore) from Bio-ORACLE and other rasters using `terra`/`sf`, and fetches species attributes from WoRMS via `worrms`. Exports `Red_sea_eDNA_enviro_var.csv`.

### Step 3 — `Red_sea_traits_clean_and_fill.Rmd` (not in the repository)
Cleans the curated trait sheet and fills missing traits from WoRMS and from relatives (≥ 80% agreement). Its output, `Red_sea_eDNA_trait_data.csv`, is committed.

### Step 4 — `Red_sea_analysis_plots_euk.Rmd`
Main analysis script. Uses non-rarefied data, combines the two primers, and counts a taxon as present in a sample at ≥ 0.1% of its reads; all analyses use presence/absence. Produces:
- bar plots (Figs 2–3) and map (Fig 1);
- multivariate GLM with `mvabund::manyglm` (Table 1, Fig 4);
- NMDS (Jaccard) with `ordisurf` GAMs (Table 2, Fig 5);
- fourth-corner and RLQ trait analysis with `ade4` (Fig 6, RLQ biplot);
- within-habitat distance decay against through-water distance (Fig 7);
- a primer / presence-cutoff sensitivity table and Venn diagrams.

### `manyglm_shrink_parallel.R`
Runs the Table 1 and sensitivity `manyglm` tests with `cor.type = "shrink"`, splitting the bootstraps across cores. Too slow for RStudio, so run it from the terminal. Reads `manyglm_inputs.rds` (written by the `export_mva_inputs` chunk) and writes `manyglm_shrink_results.rds`:

```
caffeinate -ims Rscript manyglm_shrink_parallel.R [n_cores] [n_boot] [datasets]
```

### `worms_name_check.R`
Standalone script that validates species names against WoRMS accepted names in batches. Requires internet. Run from the R console after step 1.

---

## Data files

| File | Description |
|------|-------------|
| `RS_surface_sample_data_all.xlsx` | Sample metadata (coordinates, habitat, season, depth); one sheet per primer plus a `sample_data` sheet |
| `RS_surface_filtered_data_all.csv.gz` | Filtered dataset from step 1 (reads per taxon × sample × primer), gzipped; main input for steps 2 and 4 (`data.table::fread` reads it directly; needs `R.utils`) |
| `RS_surface_dropped_samples.csv` | Every metadata sample with its status (retained / dropped and why) and reads per primer |
| `Red_sea_non_marine_taxa.csv` | Rules (rank, name, reason, source) for removing non-marine and terrestrial taxa in step 1 |
| `Red_sea_eDNA_enviro_var.csv` | Environmental covariates per sample from step 2 |
| `Red_sea_eDNA_trait_data.csv` | Trait data for all analysis taxa: taxonomy, five traits (size class, trophic level, mobility, zone, planktonic larva), WoRMS match and environment, per-trait source (literature / WoRMS / relatives / expert) and references |
| `through_water_dist_km.csv` | Pairwise through-water distances (km) between sites, from NOAA bathymetry (`marmap`) with land as a barrier; coastal sites on land cells snapped to the nearest water cell |

---

## Requirements

- R ≥ 4.3 (last run with R 4.5.3)
- Raw inputs not in this repository: DADA2 outputs (`Dada2_output/`, needed for step 1) and environmental rasters (step 2). The committed files are enough to rerun step 4 and `manyglm_shrink_parallel.R`.
- Key packages: `tidyverse`, `DECIPHER`, `decontam`, `vegan`, `ade4`, `mvabund`, `mgcv`, `terra`, `sf`, `marmap`, `worrms`, `ggplot2`, `patchwork`, `ggvenn`, `wesanderson`, `flextable`, `R.utils`

---

## Reference

Geraldi et al. (in prep.) — Red Sea coastal sediment eDNA metabarcoding.

For questions contact: nategeraldi@gmail.com
