# ============================================
# Trim-and-fill and plotting (v4 - study-level aggregation + observation-level sensitivity)
# Scheme A: inverse-variance weighted aggregation to study level by study_id, then Egger/FSN/trimfill
# Observation-level results retained as sensitivity analysis; two sets output in parallel, no overwrite
#
# Changes from v3:
# [A-1] Added aggregate_to_study(): inverse-variance weighted merge of lnRR and variance by study_id
# [A-2] Added run_bias_pipeline(): encapsulates "data -> bias -> plot -> table" into reusable pipeline,
#       called separately with level = "study" and level = "observation"
# [A-3] Output directories split into two sets: results/ (observation-level, retained) + results_study_level/ (study-level, new)
# [A-4] FSN adds optional "study-level" note; trimfill fixes side parameter (fixed in v3)
# [A-5] Row name cleanup rownames(bias_results) <- NULL
# ============================================

library(tidyverse)
library(metafor)
library(ggplot2)
library(patchwork)
library(gridExtra)

# ========== Configuration parameters ==========
DATA_DIR <- "data/raw/"

# Two output directories: observation-level (retain old results) + study-level (new)
RESULTS_DIR_OBS   <- "results/"                    # Observation-level: reuse, no overwrite
RESULTS_DIR_STUDY <- "results_study_level/"        # Study-level: new

# ========== Response variable configuration ==========
response_config <- tibble(
  file = c(
    "biomass_shoot.csv", "biomass_root.csv", "biomass_fruit.csv",
    "biomass_leave.csv", "plant_height.csv",
    "zn_shoot.csv", "zn_root.csv", "zn_plant.csv",
    "bcf_shoot.csv", "bcf_root.csv", "bcf_plant.csv",
    "part_shoot_P.csv", "part_root_P.csv"
  ),
  response_name = c(
    "Biomass_shoot", "Biomass_root", "Biomass_fruit", "Biomass_leave",
    "Plant_height", "Zn_shoot", "Zn_root", "Zn_plant",
    "BCF_shoot", "BCF_root", "BCF_plant", "P_shoot", "P_root"
  ),
  category = c(
    rep("Biomass", 5), rep("Zn_concentration", 3),
    rep("BCF", 3), rep("Phosphorus", 2)
  ),
  y_label = c(
    "Shoot Biomass (g)", "Root Biomass (g)", "Fruit Biomass (g)",
    "Leaf Biomass (g)", "Plant Height (cm)",
    "Shoot Zn (mg/kg)", "Root Zn (mg/kg)", "Whole Plant Zn (mg/kg)",
    "BCF Shoot", "BCF Root", "BCF Plant", "Shoot P (mg/g)", "Root P (mg/g)"
  )
)

# ========== Compute effect sizes (observation-level) ==========
calculate_lnRR <- function(data) {
  data %>%
    mutate(
      lnRR = log(Treatment_mean / Control_mean),
      variance = (Treatment_sd^2 / (Treatment_n * Treatment_mean^2)) +
        (Control_sd^2 / (Control_n * Control_mean^2)),
      se = sqrt(variance)
    ) %>%
    filter(is.finite(lnRR), is.finite(variance), variance > 0, se > 0,
           Treatment_mean > 0, Control_mean > 0)
}

# ========== [A-1] Inverse-variance weighted aggregation to study level by study_id ==========
# Multiple observations under same study_id:
#   w_i = 1 / variance_i
#   lnRR_study = sum(w_i * lnRR_i) / sum(w_i)
#   variance_study = 1 / sum(w_i)        (fixed-effect combined variance, standard practice)
# Also retain k_obs (number of observations contributed by study) and total_n (total sample size, for sensitivity analysis)
aggregate_to_study <- function(dat) {
  dat %>%
    mutate(w = 1 / variance) %>%
    group_by(study_id) %>%
    summarise(
      lnRR      = sum(w * lnRR) / sum(w),
      variance  = 1 / sum(w),
      k_obs     = n(),
      total_n   = sum(Treatment_n + Control_n, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    mutate(se = sqrt(variance)) %>%
    filter(is.finite(lnRR), is.finite(variance), variance > 0, se > 0)
}

# ========== Rosenthal's Fail-safe N ==========
# Note: when called on "study-level" data, k = number of studies, FSN magnitude reasonable;
#       when called on "observation-level" data, k = number of observations, FSN inflated due to repeated counting,
#       only for relative comparison, do not interpret absolute values.
calculate_fail_safe_n <- function(yi, vi, alpha = 0.05) {
  k <- length(yi)
  if (k < 2) return(NA_real_)
  z_values <- yi / sqrt(vi)
  z_combined <- sum(z_values) / sqrt(k)
  if (z_combined <= 0) return(NA_real_)
  z_alpha <- qnorm(1 - alpha)
  n_fs <- (sum(z_values)^2 / z_alpha^2) - k
  ceiling(max(0, n_fs))
}

# ========== Load observation-level data ==========
load_all_data <- function() {
  all_data <- list()
  cat("\n========== Loading observation-level data ==========\n")
  for (i in 1:nrow(response_config)) {
    file_path <- file.path(DATA_DIR, response_config$file[i])
    if (!file.exists(file_path)) {
      cat(sprintf("✗ %s: file does not exist\n", response_config$file[i])); next
    }
    dat <- tryCatch(read.csv(file_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(dat)) {
      cat(sprintf("✗ %s: read failed\n", response_config$file[i])); next
    }
    required_cols <- c("Treatment_mean", "Control_mean", "Treatment_sd",
                       "Control_sd", "Treatment_n", "Control_n", "study_id")
    missing_cols <- required_cols[!required_cols %in% colnames(dat)]
    if (length(missing_cols) > 0) {
      cat(sprintf("✗ %s: missing columns %s\n", response_config$file[i],
                  paste(missing_cols, collapse = ", "))); next
    }
    dat <- calculate_lnRR(dat)
    if (nrow(dat) >= 3) {
      all_data[[response_config$response_name[i]]] <- dat
      cat(sprintf("✓ %s: %d observations, %d studies\n",
                  response_config$response_name[i], nrow(dat),
                  length(unique(dat$study_id))))
    } else {
      cat(sprintf("✗ %s: insufficient observations (%d)\n", response_config$response_name[i], nrow(dat)))
    }
  }
  cat(sprintf("\nSuccessfully loaded %d/%d response variables\n", length(all_data), nrow(response_config)))
  all_data
}

# ========== Identify publication bias (observation-level / study-level common) ==========
identify_publication_bias <- function(processed_data, level = c("observation", "study")) {
  level <- match.arg(level)
  bias_results <- data.frame()
  
  for (resp in names(processed_data)) {
    dat <- processed_data[[resp]]
    if (is.null(dat) || nrow(dat) < 5) next
    
    # ---- Egger's test: primary standard OLS; weighted version as sensitivity ----
    egger_ols <- tryCatch({
      fit <- lm(lnRR ~ se, data = dat)
      list(p = coef(summary(fit))[2, 4], b = coef(fit)[2])
    }, error = function(e) list(p = NA_real_, b = NA_real_))
    
    egger_w <- tryCatch({
      m <- rma(yi = lnRR, vi = variance, data = dat, method = "REML")
      regtest(m, model = "rma")$pval
    }, error = function(e) NA_real_)
    
    egger_p <- egger_ols$p
    
    # ---- Fail-safe N ----
    fsn <- calculate_fail_safe_n(dat$lnRR, dat$variance, alpha = 0.05)
    n_units <- nrow(dat)                # Observation-level = observations; study-level = studies
    threshold <- 5 * n_units + 10
    fsn_robust <- if (is.na(fsn)) NA else (fsn > threshold)
    
    # ---- Trim-and-fill ----
    trimfill_result <- tryCatch({
      m <- rma(yi = lnRR, vi = variance, data = dat, method = "REML")
      tf <- trimfill(m, side = "right")   # v3 fix: side only accepts left/right
      list(
        missing_studies = tf$k0,
        adjusted_effect = exp(tf$b[1]) - 1,
        original_effect = exp(m$b[1]) - 1,
        original_p = m$pval, adjusted_p = tf$pval,
        original_ci_lb = m$ci.lb, original_ci_ub = m$ci.ub,
        adjusted_ci_lb = tf$ci.lb, adjusted_ci_ub = tf$ci.ub
      )
    }, error = function(e) list(
      missing_studies = NA_integer_, adjusted_effect = NA_real_,
      original_effect = NA_real_, original_p = NA_real_, adjusted_p = NA_real_,
      original_ci_lb = NA_real_, original_ci_ub = NA_real_,
      adjusted_ci_lb = NA_real_, adjusted_ci_ub = NA_real_
    ))
    
    # ---- Flagging logic (labels only, not formal joint test) ----
    significant_bias_egger <- !is.na(egger_p) && egger_p < 0.05
    fsn_says_unrobust  <- !is.na(fsn_robust) && !fsn_robust
    trimfill_says_asym <- !is.na(trimfill_result$missing_studies) &&
      trimfill_result$missing_studies > 0
    # Scheme 1: Egger significant AND trimfill has imputed studies -> flag
    # FSN systematically inflated due to high data precision, not used in judgment, descriptive only
    significant_bias_combined <- significant_bias_egger && trimfill_says_asym
    
    cat_idx <- which(response_config$response_name == resp)[1]
    category <- if (is.na(cat_idx)) "Unknown" else as.character(response_config$category[cat_idx])
    
    # Number of studies (observation-level data has study_id; study-level data each row one study_id)
    n_study <- if ("study_id" %in% names(dat)) length(unique(dat$study_id)) else NA_integer_
    
    bias_results <- rbind(bias_results, data.frame(
      Response = resp,
      Category = category,
      Level = level,
      N = nrow(dat),
      N_studies = n_study,
      Egger_OLS_b = round(egger_ols$b, 4),
      Egger_p = egger_p,
      Egger_weighted_p = round(egger_w, 4),
      Significant_bias_Egger = significant_bias_egger,
      Fail_safe_N = fsn,
      Threshold_5N10 = threshold,
      Exceeds_threshold = fsn_robust,
      Significant_bias_combined = significant_bias_combined,
      Missing_studies = trimfill_result$missing_studies,
      Original_effect = round(trimfill_result$original_effect * 100, 1),
      Original_CI_lower = round((exp(trimfill_result$original_ci_lb) - 1) * 100, 1),
      Original_CI_upper = round((exp(trimfill_result$original_ci_ub) - 1) * 100, 1),
      Adjusted_effect = round(trimfill_result$adjusted_effect * 100, 1),
      Adjusted_CI_lower = round((exp(trimfill_result$adjusted_ci_lb) - 1) * 100, 1),
      Adjusted_CI_upper = round((exp(trimfill_result$adjusted_ci_ub) - 1) * 100, 1),
      Original_p = trimfill_result$original_p,
      Adjusted_p = trimfill_result$adjusted_p,
      Effect_change = round((trimfill_result$adjusted_effect - trimfill_result$original_effect) * 100, 1)
    ))
  }
  rownames(bias_results) <- NULL
  bias_results
}

# ========== Funnel plot (observation-level / study-level common) ==========
plot_funnel <- function(data, response_name, category, figures_dir) {
  tryCatch({
    model <- rma(yi = lnRR, vi = variance, data = data, method = "REML")
    color_map <- c("Biomass" = "#2E86AB", "Zn_concentration" = "#A23B72",
                   "BCF" = "#F18F01", "Phosphorus" = "#4CAF50")
    category <- as.character(category)[1]
    col_single <- if (category %in% names(color_map)) unname(color_map[category]) else "#999999"
    col_vec <- rep(col_single, nrow(data))
    pch_vec <- rep(16, nrow(data))
    
    png(file.path(figures_dir, "funnel_plots", paste0(response_name, "_funnel.png")),
        width = 8, height = 6, units = "in", res = 350, bg = "white")
    funnel(model, main = paste(response_name, "Funnel Plot"),
           xlab = "Effect Size (lnRR)", ylab = "Standard Error",
           col = col_vec, pch = pch_vec, back = "white")
    dev.off()
    
    pdf(file.path(figures_dir, "funnel_plots", paste0(response_name, "_funnel.pdf")),
        width = 8, height = 6)
    funnel(model, main = paste(response_name, "Funnel Plot"),
           xlab = "Effect Size (lnRR)", ylab = "Standard Error",
           col = col_vec, pch = pch_vec)
    dev.off()
  }, error = function(e) {
    cat(sprintf("  Warning: %s funnel plot failed - %s\n", response_name, e$message))
    while (dev.cur() > 1) dev.off()
  })
}

# ========== Trim-and-fill correction comparison plot ==========
perform_trimfill_correction <- function(response_name, data, figures_dir) {
  tryCatch({
    original_model <- rma(yi = lnRR, vi = variance, data = data, method = "REML")
    trimfill_model <- trimfill(original_model, side = "right")
    results <- data.frame(
      Method = c("Original", "Trim-and-fill"),
      Effect_lnRR = round(c(original_model$b[1], trimfill_model$b[1]), 4),
      Effect_pct = round((exp(c(original_model$b[1], trimfill_model$b[1])) - 1) * 100, 1),
      CI_lower = round((exp(c(original_model$ci.lb, trimfill_model$ci.lb)) - 1) * 100, 1),
      CI_upper = round((exp(c(original_model$ci.ub, trimfill_model$ci.ub)) - 1) * 100, 1),
      P_value = round(c(original_model$pval, trimfill_model$pval), 4),
      Missing_studies = c(0, trimfill_model$k0)
    )
    cat("\n==========", response_name, "==========\n"); print(results)
    png(file.path(figures_dir, "bias_correction",
                  paste0(response_name, "_trimfill_compare.png")),
        width = 12, height = 6, units = "in", res = 350, bg = "white")
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    funnel(original_model, main = paste(response_name, "- Original"))
    funnel(trimfill_model, main = paste(response_name, "- Trim-and-fill"))
    dev.off()
    pdf(file.path(figures_dir, "bias_correction",
                  paste0(response_name, "_trimfill_compare.pdf")), width = 12, height = 6)
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    funnel(original_model, main = paste(response_name, "- Original"))
    funnel(trimfill_model, main = paste(response_name, "- Trim-and-fill"))
    dev.off()
    results
  }, error = function(e) {
    cat(sprintf("  Warning: %s trim-and-fill correction failed - %s\n", response_name, e$message))
    while (dev.cur() > 1) dev.off(); NULL
  })
}

# ========== Pre/post correction comparison forest plot ==========
plot_bias_corrected_forest <- function(bias_results, figures_dir, level_label) {
  biased_vars <- bias_results %>% filter(Significant_bias_combined == TRUE)
  if (nrow(biased_vars) == 0) { cat("\nNo variables flagged as potential bias\n"); return(NULL) }
  
  plot_data_long <- biased_vars %>%
    pivot_longer(cols = c(Original_effect, Adjusted_effect),
                 names_to = "Type", values_to = "Effect") %>%
    mutate(CI_lower = ifelse(Type == "Original_effect", Original_CI_lower, Adjusted_CI_lower),
           CI_upper = ifelse(Type == "Original_effect", Original_CI_upper, Adjusted_CI_upper),
           Method = ifelse(Type == "Original_effect", "Original", "Trim-and-fill"),
           Response_plot = factor(Response, levels = rev(unique(Response))))
  
  display_names <- c(
    "Biomass_shoot" = "Shoot Biomass", "Biomass_root" = "Root Biomass",
    "Biomass_fruit" = "Fruit Biomass", "Biomass_leave" = "Leaf Biomass",
    "Plant_height" = "Plant Height", "Zn_shoot" = "Shoot Zn",
    "Zn_root" = "Root Zn", "Zn_plant" = "Whole Plant Zn",
    "BCF_shoot" = "BCF Shoot", "BCF_root" = "BCF Root", "BCF_plant" = "BCF Plant",
    "P_shoot" = "Shoot P", "P_root" = "Root P")
  plot_data_long$Display_Name <- display_names[as.character(plot_data_long$Response)]
  plot_data_long$Response_plot <- factor(plot_data_long$Display_Name,
                                         levels = rev(unique(plot_data_long$Display_Name)))
  
  p <- ggplot(plot_data_long, aes(x = Effect, y = Response_plot, color = Method)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), width = 0.2, linewidth = 1,
                  position = position_dodge(width = 0.5)) +
    geom_point(aes(shape = Method), size = 3, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("Original" = "#2E86AB", "Trim-and-fill" = "#A23B72")) +
    scale_shape_manual(values = c("Original" = 16, "Trim-and-fill" = 17)) +
    labs(x = "Effect Size (%) [95% CI]", y = "",
         title = paste("Exploratory Trim-and-Fill Comparison (", level_label, ")", sep = ""),
         subtitle = paste("Flagged variables (Egger p < 0.05 & FSN <= 5N+10 or FSN undefined);",
                          "exploratory only, not a formal joint test")) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, color = "gray50", size = 9),
          legend.position = "bottom",
          axis.title.x = element_text(face = "bold"),
          axis.text.y = element_text(size = 10))
  
  ggsave(file.path(figures_dir, "bias_correction", "Bias_Correction_Forest.png"),
         p, width = 10, height = 6, dpi = 350, bg = "white")
  ggsave(file.path(figures_dir, "bias_correction", "Bias_Correction_Forest.pdf"),
         p, width = 10, height = 6)
  print(p); p
}

# ========== Sensitivity analysis: remove small-sample observations/studies ==========
sensitivity_remove_small_observations <- function(response_name, data, min_n = 10) {
  # Compatible with observation-level (has Treatment_n/Control_n) and study-level (only total_n)
  if (!"total_n" %in% names(data)) {
    if (all(c("Treatment_n", "Control_n") %in% names(data))) {
      data <- data %>% mutate(total_n = Treatment_n + Control_n)
    } else {
      cat(sprintf("\n%s: missing total_n or Treatment_n/Control_n, skipping sensitivity analysis\n", response_name))
      return(NULL)
    }
  }
  
  original <- rma(yi = lnRR, vi = variance, data = data, method = "REML")
  data_filtered <- data %>% filter(total_n >= min_n)
  if (nrow(data_filtered) < 3) {
    cat(sprintf("\n%s: insufficient units after removing small samples, skipping\n", response_name)); return(NULL)
  }
  filtered <- rma(yi = lnRR, vi = variance, data = data_filtered, method = "REML")
  
  has_sid <- "study_id" %in% names(data)
  n_study_orig <- if (has_sid) length(unique(data$study_id)) else NA_integer_
  n_study_filt <- if (has_sid) length(unique(data_filtered$study_id)) else NA_integer_
  
  data.frame(
    Analysis = c("Original (all units)",
                 paste0("Remove units with total_n < ", min_n)),
    N_studies = c(n_study_orig, n_study_filt),
    N_units = c(nrow(data), nrow(data_filtered)),
    Effect_pct = round(c((exp(original$b[1]) - 1) * 100,
                         (exp(filtered$b[1]) - 1) * 100), 1),
    CI_lower = round(c((exp(original$ci.lb) - 1) * 100,
                       (exp(filtered$ci.lb) - 1) * 100), 1),
    CI_upper = round(c((exp(original$ci.ub) - 1) * 100,
                       (exp(filtered$ci.ub) - 1) * 100), 1),
    P_value = round(c(original$pval, filtered$pval), 4)
  )
}

# ========== Sensitivity forest plot ==========
plot_sensitivity_forest <- function(bias_results, processed_data, figures_dir, level_label, min_n = 10) {
  biased_vars <- bias_results %>% filter(Significant_bias_combined == TRUE)
  if (nrow(biased_vars) == 0) return(NULL)
  
  sensitivity_data <- data.frame()
  for (resp in biased_vars$Response) {
    dat <- processed_data[[resp]]; if (is.null(dat)) next
    if (!"total_n" %in% names(dat)) {
      if (all(c("Treatment_n", "Control_n") %in% names(dat))) {
        dat <- dat %>% mutate(total_n = Treatment_n + Control_n)
      } else {
        next
      }
    }
    original <- rma(yi = lnRR, vi = variance, data = dat, method = "REML")
    dat_filtered <- dat %>% filter(total_n >= min_n)
    if (nrow(dat_filtered) >= 3) {
      filtered <- rma(yi = lnRR, vi = variance, data = dat_filtered, method = "REML")
      sensitivity_data <- rbind(sensitivity_data, data.frame(
        Response = resp, Type = c("Original", paste("Remove n <", min_n)),
        Effect = c((exp(original$b[1]) - 1) * 100, (exp(filtered$b[1]) - 1) * 100),
        CI_lower = c((exp(original$ci.lb) - 1) * 100, (exp(filtered$ci.lb) - 1) * 100),
        CI_upper = c((exp(original$ci.ub) - 1) * 100, (exp(filtered$ci.ub) - 1) * 100)))
    } else {
      sensitivity_data <- rbind(sensitivity_data, data.frame(
        Response = resp, Type = "Original",
        Effect = (exp(original$b[1]) - 1) * 100,
        CI_lower = (exp(original$ci.lb) - 1) * 100,
        CI_upper = (exp(original$ci.ub) - 1) * 100))
    }
  }
  if (nrow(sensitivity_data) == 0) return(NULL)
  
  display_names <- c(
    "Biomass_shoot" = "Shoot Biomass", "Biomass_root" = "Root Biomass",
    "Biomass_fruit" = "Fruit Biomass", "Biomass_leave" = "Leaf Biomass",
    "Plant_height" = "Plant Height", "Zn_shoot" = "Shoot Zn",
    "Zn_root" = "Root Zn", "Zn_plant" = "Whole Plant Zn",
    "BCF_shoot" = "BCF Shoot", "BCF_root" = "BCF Root", "BCF_plant" = "BCF Plant",
    "P_shoot" = "Shoot P", "P_root" = "Root P")
  sensitivity_data$Display_Name <- display_names[as.character(sensitivity_data$Response)]
  sensitivity_data$Response_plot <- factor(sensitivity_data$Display_Name,
                                           levels = rev(unique(sensitivity_data$Display_Name)))
  remove_label <- paste("Remove n <", min_n)
  
  p <- ggplot(sensitivity_data, aes(x = Effect, y = Response_plot, color = Type, shape = Type)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), width = 0.2, linewidth = 0.8,
                  position = position_dodge(width = 0.5)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    scale_color_manual(name = "Analysis",
                       values = c("Original" = "#2E86AB", setNames("#F18F01", remove_label))) +
    scale_shape_manual(name = "Analysis",
                       values = c("Original" = 16, setNames(17, remove_label))) +
    labs(x = "Effect Size (%) [95% CI]", y = "",
         title = paste("Sensitivity Analysis: Removing Small-Sample Units (", level_label, ")", sep = ""),
         subtitle = paste("Original effect vs. excluding units with Treatment_n + Control_n <", min_n)) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, color = "gray50", size = 10),
          legend.position = "bottom",
          axis.title.x = element_text(face = "bold"),
          axis.text.y = element_text(size = 10))
  
  ggsave(file.path(figures_dir, "bias_correction", "Sensitivity_Forest.png"),
         p, width = 10, height = 5, dpi = 350, bg = "white")
  ggsave(file.path(figures_dir, "bias_correction", "Sensitivity_Forest.pdf"),
         p, width = 10, height = 5)
  print(p); p
}

# ========== [A-2] Reusable pipeline: same logic runs observation-level / study-level ==========
run_bias_pipeline <- function(processed_data, level = c("observation", "study"),
                              results_dir, verbose = TRUE) {
  level <- match.arg(level)
  figures_dir <- file.path(results_dir, "figures")
  tables_dir  <- file.path(results_dir, "tables")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(figures_dir, "bias_correction"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(figures_dir, "funnel_plots"),    recursive = TRUE, showWarnings = FALSE)
  
  level_label <- if (level == "study") "study-level" else "observation-level"
  
  if (verbose) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════════╗\n")
    cat(sprintf("║   Publication bias assessment and trim-and-fill correction —— %-28s║\n", level_label))
    cat("╚══════════════════════════════════════════════════════════════╝\n")
  }
  
  bias_results <- identify_publication_bias(processed_data, level = level)
  write.csv(bias_results, file.path(tables_dir, "publication_bias_results.csv"),
            row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("\n========== Bias test results (%s) ==========\n", level_label))
    cols_to_print <- c("Response", "Category", "N", "N_studies",
                       "Egger_p", "Egger_weighted_p",
                       "Fail_safe_N", "Threshold_5N10", "Exceeds_threshold",
                       "Significant_bias_Egger", "Significant_bias_combined",
                       "Missing_studies")
    cols_exist <- cols_to_print[cols_to_print %in% colnames(bias_results)]
    print(bias_results[, cols_exist, drop = FALSE])
  }
  
  if (verbose) cat(sprintf("\nPlotting funnel plots (%s, PNG + PDF)...\n", level_label))
  for (resp in names(processed_data)) {
    cat_idx <- which(response_config$response_name == resp)[1]
    category <- if (is.na(cat_idx)) "Unknown" else as.character(response_config$category[cat_idx])
    plot_funnel(processed_data[[resp]], resp, category, figures_dir)
  }
  
  biased_vars <- bias_results %>% filter(Significant_bias_combined == TRUE)
  
  if (nrow(biased_vars) > 0) {
    if (verbose) cat(sprintf("\n========== Running trim-and-fill analysis for flagged variables (%s) ==========\n", level_label))
    trimfill_results <- list()
    for (i in 1:nrow(biased_vars)) {
      resp <- biased_vars$Response[i]
      trimfill_results[[resp]] <- perform_trimfill_correction(resp, processed_data[[resp]], figures_dir)
    }
    
    if (verbose) cat(sprintf("\nPlotting pre/post correction comparison forest plot (%s)...\n", level_label))
    plot_bias_corrected_forest(bias_results, figures_dir, level_label)
    
    if (verbose) cat(sprintf("\nPlotting sensitivity analysis forest plot (%s)...\n", level_label))
    plot_sensitivity_forest(bias_results, processed_data, figures_dir, level_label, min_n = 10)
    
    if (verbose) cat(sprintf("\n========== Sensitivity analysis: removing observations/studies with total sample size <10 (%s) ==========\n", level_label))
    sensitivity_results <- list()
    for (i in 1:nrow(biased_vars)) {
      resp <- biased_vars$Response[i]
      sensitivity_results[[resp]] <- sensitivity_remove_small_observations(
        resp, processed_data[[resp]], min_n = 10)
      if (!is.null(sensitivity_results[[resp]]) && verbose) print(sensitivity_results[[resp]])
    }
    sensitivity_df <- bind_rows(lapply(names(sensitivity_results), function(x) {
      if (!is.null(sensitivity_results[[x]])) {
        df <- sensitivity_results[[x]]; df$Response <- x; return(df)
      }
      return(NULL)
    }))
    if (!is.null(sensitivity_df) && nrow(sensitivity_df) > 0) {
      write.csv(sensitivity_df,
                file.path(tables_dir, "sensitivity_remove_small_units.csv"),
                row.names = FALSE)
    }
  } else {
    if (verbose) cat(sprintf("\n✅ No variables flagged as potential bias (%s)\n", level_label))
    egger_biased <- bias_results %>% filter(Significant_bias_Egger == TRUE)
    if (nrow(egger_biased) > 0 && verbose) {
      cat("\nNote: following variables have significant Egger but FSN > threshold, effect robust:\n")
      print(egger_biased[, c("Response", "Egger_p", "Fail_safe_N", "Threshold_5N10")])
    }
  }
  
  if (verbose) {
    cat("\n========================================\n")
    cat(sprintf("Complete (%s)! Output directory: %s\n", level_label, results_dir))
    cat("========================================\n")
  }
  
  list(data = processed_data, bias_results = bias_results, biased_vars = biased_vars)
}

# ============================================
# Main workflow
# ============================================
cat("\n══════════════════════════════════════════════════════════════╗\n")
cat("║   Publication bias assessment v4 (study-level primary + observation-level sensitivity, two sets in parallel)     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")

# 1. Load observation-level data
processed_obs <- load_all_data()
if (length(processed_obs) == 0) stop("No data successfully loaded!")

# 2. Aggregate to study level
cat("\n========== Inverse-variance weighted aggregation to study level by study_id ==========\n")
processed_study <- list()
for (resp in names(processed_obs)) {
  dat_study <- aggregate_to_study(processed_obs[[resp]])
  if (nrow(dat_study) >= 5) {   # Egger needs >=5 units to be meaningful
    processed_study[[resp]] <- dat_study
    cat(sprintf("✓ %s: %d observations -> %d studies\n",
                resp, nrow(processed_obs[[resp]]), nrow(dat_study)))
  } else {
    cat(sprintf("✗ %s: insufficient studies after aggregation (%d), skipping study-level analysis\n", resp, nrow(dat_study)))
  }
}

# 3. Study-level analysis (primary results)
res_study <- run_bias_pipeline(processed_study, level = "study",
                               results_dir = RESULTS_DIR_STUDY, verbose = TRUE)

# 4. Observation-level analysis (sensitivity / appendix; reuse old directory, does not overwrite old files but refreshes)
res_obs <- run_bias_pipeline(processed_obs, level = "observation",
                             results_dir = RESULTS_DIR_OBS, verbose = TRUE)

# 5. Side-by-side comparison table (for direct citation in paper)
cat("\n========== Study-level vs observation-level comparison ==========\n")
compare <- full_join(
  res_study$bias_results %>%
    select(Response, N_study = N, N_studies_study = N_studies,
           Egger_p_study = Egger_p, FSN_study = Fail_safe_N,
           Thr_study = Threshold_5N10, Miss_study = Missing_studies,
           Combined_study = Significant_bias_combined),
  res_obs$bias_results %>%
    select(Response, N_obs = N, N_studies_obs = N_studies,
           Egger_p_obs = Egger_p, FSN_obs = Fail_safe_N,
           Thr_obs = Threshold_5N10, Miss_obs = Missing_studies,
           Combined_obs = Significant_bias_combined),
  by = "Response"
) %>% arrange(Response)
print(compare)
write.csv(compare, file.path(RESULTS_DIR_STUDY, "tables",
                             "study_vs_observation_comparison.csv"),
          row.names = FALSE, na = "")

cat("\n========================================\n")
cat("All complete! Two sets of results output in parallel:\n")
cat("  - Study-level (primary results):", RESULTS_DIR_STUDY, "\n")
cat("  - Observation-level (sensitivity/appendix):", RESULTS_DIR_OBS, "\n")
cat("  - Comparison table:", file.path(RESULTS_DIR_STUDY, "tables", "study_vs_observation_comparison.csv"), "\n")
cat("========================================\n")

# ============================================
# Additional feature v2: redraw all funnel plots and combine into 4×4 panel
# - Maintain original 8:6 aspect ratio, each cell 4:3
# - No additional titles (original figure already has main)
# - Output PNG (one per page) + PDF (multi-page merged into one file)
# ============================================

combine_funnel_plots_v2 <- function(processed_data, figures_dir,
                                    page_size = 16, ncol = 4,
                                    cell_w = 4, cell_h = 3,   # Inches per cell, 4:3
                                    out_name = "Funnel_Plots_Combined") {
  funnel_dir <- file.path(figures_dir, "funnel_plots")
  
  # Keep only response variables with actual data, in order consistent with response_config
  resps <- response_config$response_name
  resps <- resps[resps %in% names(processed_data)]
  if (length(resps) == 0) {
    cat("⚠ No available data, skipping combination\n"); return(invisible(NULL))
  }
  
  n_page <- ceiling(length(resps) / page_size)
  nrow_pg <- ceiling(page_size / ncol)
  page_w <- cell_w * ncol
  page_h <- cell_h * nrow_pg
  
  # Plotting function: draw i-th funnel plot on current device
  draw_one <- function(resp) {
    dat <- processed_data[[resp]]
    model <- tryCatch(
      rma(yi = lnRR, vi = variance, data = dat, method = "REML"),
      error = function(e) NULL)
    if (is.null(model)) { plot.new(); return(invisible()) }
    
    cat_idx <- which(response_config$response_name == resp)[1]
    category <- if (is.na(cat_idx)) "Unknown" else
      as.character(response_config$category[cat_idx])
    color_map <- c("Biomass" = "#2E86AB", "Zn_concentration" = "#A23B72",
                   "BCF" = "#F18F01", "Phosphorus" = "#4CAF50")
    col_single <- if (category %in% names(color_map))
      unname(color_map[category]) else "#999999"
    
    funnel(model,
           main = paste(resp, "Funnel Plot"),
           xlab = "Effect Size (lnRR)", ylab = "Standard Error",
           col = rep(col_single, nrow(dat)), pch = rep(16, nrow(dat)),
           back = "gray95",cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.7)
  }
  
  # ---- Output PNG: one per page ----
  png_files <- character(0)
  for (pg in seq_len(n_page)) {
    idx <- ((pg - 1) * page_size + 1):min(pg * page_size, length(resps))
    idx <- idx[!is.na(idx)]
    png_out <- file.path(funnel_dir, sprintf("%s_page%d.png", out_name, pg))
    png(png_out, width = page_w, height = page_h, units = "in",
        res = 300, bg = "white")
    par(mfrow = c(nrow_pg, ncol), mar = c(3, 3, 2, 1), oma = c(0, 0, 2, 0))
    for (i in idx) draw_one(resps[i])
    # Empty positions remain blank (mfrow automatically leaves blank, no extra handling needed)
    mtext(sprintf("Funnel Plots (page %d/%d)", pg, n_page),
          outer = TRUE, cex = 1.2, font = 2)
    dev.off()
    png_files <- c(png_files, png_out)
    cat(sprintf("  ✓ PNG page %d/%d: %s\n", pg, n_page, basename(png_out)))
  }
  
  # ---- Output PDF: multi-page merged into one file ----
  pdf_out <- file.path(funnel_dir, paste0(out_name, ".pdf"))
  pdf(pdf_out, width = page_w, height = page_h, onefile = TRUE)
  for (pg in seq_len(n_page)) {
    idx <- ((pg - 1) * page_size + 1):min(pg * page_size, length(resps))
    idx <- idx[!is.na(idx)]
    par(mfrow = c(nrow_pg, ncol), mar = c(3, 3, 2, 1), oma = c(0, 0, 2, 0))
    for (i in idx) draw_one(resps[i])
    mtext(sprintf("Funnel Plots (page %d/%d)", pg, n_page),
          outer = TRUE, cex = 1.2, font = 2)
  }
  dev.off()
  cat(sprintf("  ✓ PDF merged file: %s\n", basename(pdf_out)))
  
  invisible(c(png_files, pdf_out))
}
# 6. Combine funnel plots 4×4 (redraw, preserve ratio, PNG + multi-page PDF)
cat("\n========== Combining funnel plots (4×4) ==========\n")

combine_funnel_plots_v2(processed_study, file.path(RESULTS_DIR_STUDY, "figures"),
                        page_size = 16, ncol = 4, cell_w = 4, cell_h = 3,
                        out_name = "Funnel_Plots_Combined_study")

combine_funnel_plots_v2(processed_obs, file.path(RESULTS_DIR_OBS, "figures"),
                        page_size = 16, ncol = 4, cell_w = 4, cell_h = 3,
                        out_name = "Funnel_Plots_Combined_observation")