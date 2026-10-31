# ============================================
# Module 2: Overall effect size estimation (supports Zn_Level grouping) - fixed version
# ============================================

source("00_config.R")

# ============================================
# Function: manually compute I²
# ============================================
calc_I2_for_rma_mv <- function(model, dat = NULL) {
  if(is.null(model)) return(NA)
  
  if(!is.null(model$I2) && length(model$I2) > 0) {
    if(!is.na(model$I2) && is.finite(model$I2) && model$I2 > 0) {
      return(round(model$I2, 1))
    }
  }
  
  if(!is.null(model$QE) && !is.null(model$k) && !is.null(model$p)) {
    Q <- model$QE
    df <- model$k - model$p
    if(Q > 0 && df > 0) {
      I2 <- ((Q - df) / Q) * 100
      I2 <- max(0, I2)
      return(round(I2, 1))
    }
  }
  
  tau2 <- model$tau2
  if(!is.null(model$sigma2) && length(model$sigma2) > 0) {
    sigma2_avg <- mean(model$sigma2[is.finite(model$sigma2) & model$sigma2 > 0], na.rm = TRUE)
    if(!is.na(sigma2_avg) && sigma2_avg > 0 && is.finite(sigma2_avg)) {
      total_var <- tau2 + sigma2_avg
      if(total_var > 0 && is.finite(total_var)) {
        I2 <- 100 * tau2 / total_var
        return(round(I2, 1))
      }
    }
  }
  
  if(!is.null(dat) && nrow(dat) > 0) {
    if("variance" %in% colnames(dat) && "lnRR" %in% colnames(dat)) {
      vi <- dat$variance
      yi <- dat$lnRR
      weights <- 1 / vi
      y_bar <- sum(weights * yi) / sum(weights)
      Q_fixed <- sum(weights * (yi - y_bar)^2)
      df_fixed <- length(yi) - 1
      if(Q_fixed > 0 && df_fixed > 0) {
        I2 <- ((Q_fixed - df_fixed) / Q_fixed) * 100
        I2 <- max(0, I2)
        return(round(I2, 1))
      }
    }
  }
  
  return(0)
}

apply_variance_floor <- function(dat, floor_value = 0.001) {
  if(is.null(dat)) return(NULL)
  n_modified <- sum(dat$variance < floor_value, na.rm = TRUE)
  dat <- dat %>%
    mutate(
      variance = pmax(variance, floor_value),
      weight = 1 / variance,
      pct_lower = (exp(lnRR - 1.96 * sqrt(variance)) - 1) * 100,
      pct_upper = (exp(lnRR + 1.96 * sqrt(variance)) - 1) * 100
    )
  if(n_modified > 0) {
    cat(sprintf("    ✓ Corrected %d observations with variance < %.3f\n", n_modified, floor_value))
  }
  return(dat)
}

# ============================================
# Forest plot function (supports grouping) - fixed version
# ============================================
plot_forest_grouped <- function(dat, response_name, category, y_label, 
                                overall_pct, overall_pct_lower, overall_pct_upper,
                                zn_level = NULL, group_label = "All") {
  
  forest_data <- dat %>%
    mutate(
      study_label = paste0(study_id, " (n = ", Treatment_n + Control_n, ")"),
      pct = pct_change,
      pct_lower = pct_lower,
      pct_upper = pct_upper
    ) %>%
    arrange(desc(pct))
  
  if(nrow(forest_data) > 50) {
    forest_data <- forest_data %>% slice_head(n = 50)
    cat("  Note: forest plot shows only first 50 observations\n")
  }
  
  # Select color
  if(!is.null(zn_level) && zn_level %in% names(ZN_LEVEL_COLORS)) {
    color <- ZN_LEVEL_COLORS[zn_level]
  } else {
    color <- ifelse(category == "Biomass", "#2E86AB",
                    ifelse(category == "Zn_concentration", "#A23B72",
                           ifelse(category == "BCF", "#F18F01",
                                  ifelse(category == "Phosphorus", "#4CAF50", "gray70"))))
  }
  
  p <- ggplot(forest_data, aes(x = reorder(study_label, pct), y = pct)) +
    geom_point(size = 2.5, color = color) +
    geom_errorbar(aes(ymin = pct_lower, ymax = pct_upper), 
                  width = 0.3, color = color, linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    coord_flip() +
    labs(
      title = paste(response_name, "-", group_label),
      subtitle = paste("Overall effect:", sprintf("%.2f%% [%.2f, %.2f]", 
                                                  overall_pct, overall_pct_lower, overall_pct_upper)),
      x = "Study / Observation",
      y = paste("Change (%) -", y_label)
    ) +
    THEME +
    theme(axis.text.y = element_text(size = 7))
  
  # ========== Fix: safely create directory ==========
  if(SAVE_PLOTS) {
    # Build forest plot directory path
    if(!is.null(zn_level)) {
      forest_dir <- file.path(FIGURES_DIR, "forest_plots", zn_level)
    } else {
      forest_dir <- file.path(FIGURES_DIR, "forest_plots", "all")
    }
    
    # Safely create directory
    if(!dir.exists(forest_dir)) {
      dir.create(forest_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Save plot
    file_path <- file.path(forest_dir, paste0(response_name, "_forest.png"))
    ggsave(file_path, p, width = 10, height = max(6, nrow(forest_data) * 0.2), 
           dpi = 300, bg = "white")
    cat("  ✓ Forest plot saved:", file_path, "\n")
  }
  
  return(p)
}

# ============================================
# Main function: run overall multilevel model (supports grouping)
# ============================================
run_overall_analysis <- function(dat, response_name, category, y_label, 
                                 zn_level = NULL, group_label = "All",
                                 variance_floor = 0.001, 
                                 verbose = FALSE) {
  
  cat(sprintf("\n--- %s (%s) ---\n", response_name, group_label))
  
  if(nrow(dat) < 3) {
    cat("  Insufficient sample size (n =", nrow(dat), "), skipping\n")
    return(NULL)
  }
  
  n_obs_original <- nrow(dat)
  n_studies_original <- n_distinct(dat$study_id)
  
  cat(sprintf("  Variance range (before correction): [%.6f, %.4f]\n", 
              min(dat$variance, na.rm = TRUE), 
              max(dat$variance, na.rm = TRUE)))
  
  dat <- apply_variance_floor(dat, variance_floor)
  
  cat(sprintf("  Variance range (after correction): [%.6f, %.4f]\n", 
              min(dat$variance, na.rm = TRUE), 
              max(dat$variance, na.rm = TRUE)))
  
  # Multilevel random-effects model
  model <- tryCatch({
    rma.mv(
      yi = lnRR,
      V = variance,
      random = RANDOM_STRUCTURE,
      data = dat,
      method = "REML"
    )
  }, error = function(e) {
    cat("  Model fitting failed:", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(model)) return(NULL)
  
  if(verbose) {
    cat("\n=== Model diagnostics ===\n")
    cat("  k (observations):", model$k, "\n")
    cat("  p (parameters):", model$p, "\n")
    cat("  tau2:", model$tau2, "\n")
    cat("  sigma2:", paste(model$sigma2, collapse = ", "), "\n")
    cat("  QE:", model$QE, "\n")
    cat("  I2 (model built-in):", model$I2, "\n")
    cat("================\n")
  }
  
  # Extract results
  overall_lnRR <- model$b[1]
  overall_lower <- model$ci.lb
  overall_upper <- model$ci.ub
  overall_pct <- (exp(overall_lnRR) - 1) * 100
  overall_pct_lower <- (exp(overall_lower) - 1) * 100
  overall_pct_upper <- (exp(overall_upper) - 1) * 100
  
  q_df <- model$k - model$p
  q_ratio <- model$QE / q_df
  I2_value <- calc_I2_for_rma_mv(model, dat)
  
  cat(sprintf("  Overall effect: %.2f%% [95%% CI: %.2f%%, %.2f%%]\n",
              overall_pct, overall_pct_lower, overall_pct_upper))
  cat(sprintf("  Q = %.2f (df = %d, p = %.4f)\n", 
              model$QE, q_df, model$QEp))
  cat(sprintf("  Q/df = %.2f\n", q_ratio))
  cat(sprintf("  I² = %.1f%%\n", I2_value))
  cat(sprintf("  tau² = %.3f\n", model$tau2))
  
  if(!is.na(I2_value) && is.finite(I2_value)) {
    if(I2_value > 75) {
      cat("  ⚠ I² > 75%, high heterogeneity\n")
    } else if(I2_value > 50) {
      cat("  ⚠ I² > 50%, moderate heterogeneity\n")
    }
  }
  
  # Forest plot
  if(SAVE_PLOTS && nrow(dat) <= 100) {
    tryCatch({
      plot_forest_grouped(
        dat = dat,
        response_name = response_name,
        category = category,
        y_label = y_label,
        overall_pct = overall_pct,
        overall_pct_lower = overall_pct_lower,
        overall_pct_upper = overall_pct_upper,
        zn_level = zn_level,
        group_label = group_label
      )
    }, error = function(e) {
      cat("  Forest plot drawing failed:", e$message, "\n")
    })
  }
  
  return(list(
    response = response_name,
    category = category,
    zn_level = zn_level,
    group_label = group_label,
    model = model,
    n_obs = nrow(dat),
    n_obs_original = n_obs_original,
    n_studies = n_distinct(dat$study_id),
    n_studies_original = n_studies_original,
    overall_lnRR = overall_lnRR,
    overall_pct = overall_pct,
    overall_pct_lower = overall_pct_lower,
    overall_pct_upper = overall_pct_upper,
    I2 = I2_value,
    Q = model$QE,
    Q_df = q_df,
    Q_p = model$QEp,
    Q_ratio = q_ratio,
    variance_floor = variance_floor,
    tau2 = model$tau2,
    pval = model$pval
  ))
}

# ============================================
# Batch run all response variables (supports grouping)
# ============================================
run_all_overall <- function(processed_data, variance_floor = 0.001, 
                            zn_level = NULL, group_label = "All",
                            verbose = FALSE) {
  
  cat("\n========================================\n")
  cat(sprintf("  Overall effect analysis - %s\n", group_label))
  cat("========================================\n")
  cat(sprintf("  Variance floor: %.3f\n", variance_floor))
  cat(sprintf("  Random-effects structure: within-study + within-observation\n"))
  cat("========================================\n")
  
  results <- list()
  
  for(i in 1:nrow(RESPONSE_FILES)) {
    resp <- RESPONSE_FILES$response_name[i]
    dat <- processed_data[[resp]]
    
    if(is.null(dat)) {
      cat(sprintf("\n⚠ %s: data not found, skipping\n", resp))
      next
    }
    
    result <- run_overall_analysis(
      dat = dat, 
      response_name = resp, 
      category = RESPONSE_FILES$category[i],
      y_label = RESPONSE_FILES$y_label[i],
      zn_level = zn_level,
      group_label = group_label,
      variance_floor = variance_floor,
      verbose = verbose
    )
    
    if(!is.null(result)) {
      results[[resp]] <- result
    }
  }
  
  # Summary table
  if(length(results) > 0) {
    summary_table <- map_dfr(results, ~{
      tibble(
        Response = .x$response,
        Category = .x$category,
        Zn_Level = ifelse(is.null(.x$zn_level), "All", .x$zn_level),
        N_obs = .x$n_obs,
        N_studies = .x$n_studies,
        Effect_pct = sprintf("%.2f [%.2f, %.2f]", 
                             .x$overall_pct, .x$overall_pct_lower, .x$overall_pct_upper),
        I2 = sprintf("%.1f%%", .x$I2),
        Q = sprintf("%.2f", .x$Q),
        Q_p = sprintf("%.4f", .x$Q_p),
        tau2 = sprintf("%.3f", .x$tau2)
      )
    })
    
    if(SAVE_TABLES) {
      if(!dir.exists(TABLES_DIR)) dir.create(TABLES_DIR, recursive = TRUE)
      if(is.null(zn_level)) {
        fname <- "overall_effects_all.csv"
      } else {
        fname <- paste0("overall_effects_", zn_level, ".csv")
      }
      write_csv(summary_table, file.path(TABLES_DIR, fname))
      cat("\n✅ Summary table saved:", file.path(TABLES_DIR, fname), "\n")
    }
    
    cat("\n========================================\n")
    cat(sprintf("  Overall effect summary - %s\n", group_label))
    cat("========================================\n")
    print(summary_table)
  }
  
  return(results)
}

# ============================================
# Run main program
# ============================================
if(exists("processed_data")) {
  
  # Get Zn_Level column name (for checking)
  has_zn_level <- FALSE
  for(resp in names(processed_data)) {
    dat <- processed_data[[resp]]
    if("Zn_Level" %in% colnames(dat)) {
      has_zn_level <- TRUE
      break
    }
  }
  
  if(!has_zn_level) {
    cat("⚠️ Zn_Level column not found in data, will run overall analysis only\n")
  }
  
  # All data
  overall_results_all <- run_all_overall(
    processed_data = processed_data, 
    zn_level = NULL,
    group_label = "All",
    variance_floor = 0.001,
    verbose = FALSE
  )
  saveRDS(overall_results_all, file.path(PROCESSED_DIR, "overall_results_all.rds"))
  cat("\n✅ Overall results saved: overall_results_all.rds\n")
  
  # If grouping enabled and Zn_Level column exists
  if(USE_ZN_LEVEL_GROUP && has_zn_level) {
    
    # Try loading grouped data
    low_file <- file.path(PROCESSED_DIR, "all_processed_data_low.rds")
    high_file <- file.path(PROCESSED_DIR, "all_processed_data_high.rds")
    
    if(file.exists(low_file) && file.exists(high_file)) {
      
      # Low Zn group
      processed_data_low <- readRDS(low_file)
      if(length(processed_data_low) > 0) {
        overall_results_low <- run_all_overall(
          processed_data = processed_data_low,
          zn_level = "low",
          group_label = "Low Zn Level",
          variance_floor = 0.001,
          verbose = FALSE
        )
        saveRDS(overall_results_low, file.path(PROCESSED_DIR, "overall_results_low.rds"))
        cat("✅ Saved: overall_results_low.rds\n")
      }
      
      # High Zn group
      processed_data_high <- readRDS(high_file)
      if(length(processed_data_high) > 0) {
        overall_results_high <- run_all_overall(
          processed_data = processed_data_high,
          zn_level = "high",
          group_label = "High Zn Level",
          variance_floor = 0.001,
          verbose = FALSE
        )
        saveRDS(overall_results_high, file.path(PROCESSED_DIR, "overall_results_high.rds"))
        cat("✅ Saved: overall_results_high.rds\n")
      }
      
      cat("\n========================================\n")
      cat("  ✅ Overall effect analysis complete!\n")
      cat("========================================\n")
      cat("Output files:\n")
      cat("  - overall_results_all.rds (all data)\n")
      cat("  - overall_results_low.rds (Low Zn group)\n")
      cat("  - overall_results_high.rds (High Zn group)\n")
      cat("  - results/tables/overall_effects_all.csv\n")
      cat("  - results/tables/overall_effects_low.csv\n")
      cat("  - results/tables/overall_effects_high.csv\n")
      cat("  - results/figures/forest_plots/all/ (forest plots)\n")
      cat("  - results/figures/forest_plots/low/ (forest plots)\n")
      cat("  - results/figures/forest_plots/high/ (forest plots)\n")
      
    } else {
      cat("⚠️ Grouped data files do not exist, please run 00_winsorize_data.R first\n")
    }
  }
  
} else {
  cat("Error: processed_data does not exist, please run 01_effect_size.R or 00_winsorize_data.R first\n")
}