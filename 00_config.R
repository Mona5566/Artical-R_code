# ============================================
# Configuration file: 00_config.R
# Paths, parameters, variable lists
# ============================================

library(tidyverse)
library(metafor)
library(randomForest)
library(segmented)
library(ggplot2)
library(patchwork)

# ========== Data source selection ==========
USE_WINSORIZED_DATA <- TRUE # TRUE: use winsorized data, FALSE: use raw data

# ========== Zn_Level grouping settings ==========
# TRUE: analyze by Zn_Level (low/high) groups
# FALSE: do not group, analyze all data
USE_ZN_LEVEL_GROUP <- TRUE

# ========== Path settings ==========
if(USE_WINSORIZED_DATA) {
  DATA_DIR <- "data/winsorized/"
  cat("\n>>> Using winsorized data <<<\n")
} else {
  DATA_DIR <- "data/raw/"
  cat("\n>>> Using raw data <<<\n")
}

PROCESSED_DIR <- "data/processed/"
RESULTS_DIR <- "results/"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
TABLES_DIR <- file.path(RESULTS_DIR, "tables")

# Create directories
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, showWarnings = FALSE)
dir.create(RESULTS_DIR, showWarnings = FALSE)
dir.create(FIGURES_DIR, showWarnings = FALSE)
dir.create(TABLES_DIR, showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "forest_plots"), showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "rf_importance"), showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "threshold_plots"), showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "metareg_plots"), showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "funnel_plots"), showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "loo_plots"), showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "subgroup_plots"), showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "pdp_plots"), showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "grouped_importance"), showWarnings = FALSE)

# ========== Response variable list (13) ==========
RESPONSE_FILES <- tibble(
  file = c(
    "biomass_shoot.csv",
    "biomass_root.csv",
    "biomass_fruit.csv",
    "biomass_leave.csv",
    "plant_height.csv",
    "zn_shoot.csv",
    "zn_root.csv",
    "zn_plant.csv",
    "bcf_shoot.csv",
    "bcf_root.csv",
    "bcf_plant.csv",
    "part_shoot_P.csv",
    "part_root_P.csv"
  ),
  response_name = c(
    "Biomass_shoot",
    "Biomass_root",
    "Biomass_fruit",
    "Biomass_leave",
    "Plant_height",
    "Zn_shoot",
    "Zn_root",
    "Zn_plant",
    "BCF_shoot",
    "BCF_root",
    "BCF_plant",
    "P_shoot",
    "P_root"
  ),
  category = c(
    rep("Biomass", 5),
    rep("Zn_concentration", 3),
    rep("BCF", 3),
    rep("Phosphorus", 2)
  ),
  y_label = c(
    "Biomass (g)",
    "Biomass (g)",
    "Biomass (g)",
    "Biomass (g)",
    "Plant height (cm)",
    "Zn concentration (mg/kg)",
    "Zn concentration (mg/kg)",
    "Zn concentration (mg/kg)",
    "Bioconcentration Factor",
    "Bioconcentration Factor",
    "Bioconcentration Factor",
    "P concentration (mg/g)",
    "P concentration (mg/g)"
  )
)

# ========== Moderators ==========
CONTINUOUS_MODERATORS <- c("Total_Zn", "soil_pH")
CATEGORICAL_MODERATORS <- c(
  "Type_soil", "Type_plant", "Type_AMF2", "Soil_P", 
  "Type_root", "Type_pH", "Duration"
)
ALL_MODERATORS <- c(CONTINUOUS_MODERATORS, CATEGORICAL_MODERATORS)

# ========== Zn_Level group labels ==========
ZN_LEVEL_LABELS <- c("low" = "Low Zn Level", "high" = "High Zn Level")

# ========== Analysis parameters ==========
RANDOM_STRUCTURE <- list(~ 1 | study_id, ~ 1 | obs_id)
RF_NTREE <- 500
THRESHOLD_MIN_OBS <- 30
SIG_LEVEL <- 0.05

# ========== Winsorization parameters ==========
WINSORIZE <- TRUE
WINSORIZE_PERCENTILE <- 0.95

# ========== Plotting parameters ==========
SAVE_PLOTS <- TRUE
SAVE_TABLES <- TRUE

THEME <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5)
  )

COLOR_SCHEME <- c(
  "Biomass" = "#2E86AB", 
  "Zn_concentration" = "#A23B72", 
  "BCF" = "#F18F01",
  "Phosphorus" = "#4CAF50"
)

# Zn_Level colors
ZN_LEVEL_COLORS <- c("low" = "#2E86AB", "high" = "#A23B72")

get_response_color <- function(category) { return(COLOR_SCHEME[category]) }
get_response_label <- function(response_name) {
  label <- RESPONSE_FILES$y_label[RESPONSE_FILES$response_name == response_name]
  if(length(label) == 0) return("Value")
  return(label)
}
get_response_category <- function(response_name) {
  cat <- RESPONSE_FILES$category[RESPONSE_FILES$response_name == response_name]
  if(length(cat) == 0) return(NA)
  return(cat)
}

print_config_summary <- function() {
  cat("\n========== Configuration Summary ==========\n")
  cat("Data source:", ifelse(USE_WINSORIZED_DATA, "Winsorized data", "Raw data"), "\n")
  cat("Winsorization percentile:", WINSORIZE_PERCENTILE * 100, "%\n")
  cat("Zn_Level grouping:", ifelse(USE_ZN_LEVEL_GROUP, "Yes (low/high)", "No"), "\n")
  cat("Number of response variables:", nrow(RESPONSE_FILES), "\n")
  cat("Save plots:", SAVE_PLOTS, "\n")
  cat("Save tables:", SAVE_TABLES, "\n")
  cat("===============================\n")
}

if(interactive()) print_config_summary()