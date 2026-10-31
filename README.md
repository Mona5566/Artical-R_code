# Artical-R_code
AMF-Zn-meta
Code Documentation and Operation Manual
1. Project Overview
This is a comprehensive meta-analysis pipeline investigating the effects of arbuscular mycorrhizal fungi (AMF) on plant responses under different zinc (Zn) levels. The analysis is organized into modular R scripts that handle data preprocessing, effect size calculation, overall and subgroup analyses, meta-regression, variable importance analysis, publication bias assessment, and visualization. The pipeline supports three analysis groups: All (combined data), Low Zn, and High Zn subgroups, and covers 13 response variables across four categories: Biomass, Zn concentration, BCF (Bioconcentration Factor), and Phosphorus.
2. Directory Structure
project_root/
├── 00_config.R                  # Configuration: paths, parameters, variable lists
├── 00_winsorize_data.R          # Winsorization of lnRR effect sizes + grouping
├── 01_effect_size.R             # Effect size (lnRR) computation
├── 02_overall_analysis.R        # Overall multilevel meta-analysis
├── Fig1.R                       # Global map of study sites
├── Fig2.R                       # Low vs High Zn comparison (forest, scatter, heatmap)
├── Fig3.R                       # Three-level meta + plant part comparison
├── Fig4.R                       # Random forest variable importance
├── Fig5.1.R                     # Smart meta-regression model selection + summary table
├── Fig5.2.R                     # Meta-regression fit plots by best model
├── FigS2 ... Egger+Fail-safeN.R # Publication bias assessment (study & observation level)
├── FigS3.R                      # Meta-regression combined panels (Total_Zn, soil_pH)
├── FigS4.R                      # Subgroup analysis dot plots
├── FigS5-S7subgroup_forest.R    # Subgroup forest plots (three-level model)
├── Table1Meta.R                 # Meta-regression results tables
│
├── data/
│   ├── raw/                     # Raw CSV input (13 files + mapdata2.csv)
│   ├── winsorized/              # Winsorized CSV output + log
│   └── processed/               # RDS outputs (all_processed_data*.rds, results)
│
├── results/
│   ├── figures/                 # All plots (PNG + PDF)
│   ├── tables/                  # CSV summary tables
│   └── results_study_level/     # Study-level bias analysis outputs
│
└── output_figures/              # Figure 3/4 specific outputs


3.Statistical Methods Summary
(1)Effect size: log response ratio, lnRR = log(Treatment_mean / Control_mean), with variance from delta method.
(2)Meta-analytic models: multilevel random-effects (rma.mv) with observation nested within study (~1|study_id/obs_id); REML estimation.
(3)Heterogeneity: τ² per level, I² computed as variance-component proportion (custom function in 02_overall_analysis.R).
(4)Group comparison (Low vs High): meta-regression with mods = ~group using a two-level model; also a Welch t-test fallback.
(5)Meta-regression: automatic model selection among linear, quadratic (centered), log, and segmented models by AIC; significant models preferred.
(6)Multiple testing: Benjamini-Hochberg FDR applied to subgroup Q_M tests and pairwise contrasts.
(7)Publication bias: Egger’s regression (OLS + weighted sensitivity), Rosenthal’s fail-safe N (threshold 5k+10), trim-and-fill (side = "right"). Results reported at study level (primary) and observation level (sensitivity).
(8)Variable importance: random forest with %IncMSE, normalized within each model; group contributions averaged across responses without re-normalization.
