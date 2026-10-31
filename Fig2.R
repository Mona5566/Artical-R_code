# ============================================
# Module: Zn_Level two-group comparison analysis (earth-tone color scheme)
# Compare AMF effect differences between Low Zn vs High Zn
# Statistical method: two-level random-effects meta-regression + FDR correction + prediction intervals
# Visualization: earth-tone color scheme (Nature/Ecology style)
# ============================================

# ============================================
# 0. Environment setup
# ============================================

source("00_config.R")

# Load required packages
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(scales)

# Check and load metafor
if (!requireNamespace("metafor", quietly = TRUE)) {
  install.packages("metafor")
}
library(metafor)

# ============================================
# 0.5 Resolve function conflicts
# ============================================

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
rename <- dplyr::rename
pivot_longer <- tidyr::pivot_longer
pivot_wider <- tidyr::pivot_wider

# Key fix: margin / unit masked by patchwork etc.
# causing margin.default(t = 6, ...) to throw "missing argument observed" error
margin <- ggplot2::margin
unit <- ggplot2::unit

options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)

# ============================================
# 0.7 Statistical configuration
# ============================================

USE_IMPROVED_STATS <- TRUE
P_ADJ_METHOD <- "fdr"
CI_LEVEL <- 0.95

# ============================================
# 1. Top-journal style theme + earth-tone colors
# ============================================

theme_journal <- function(base_size = 10, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = base_size),
      axis.text.x = element_text(size = base_size, margin = margin(t = 2)),
      axis.text.y = element_text(size = base_size, margin = margin(r = 2)),
      axis.title = element_text(size = base_size + 1, face = "plain"),
      axis.title.x = element_text(margin = margin(t = 6, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.1, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = base_size - 0.5),
      legend.title = element_text(size = base_size, face = "plain"),
      legend.position = "top",
      legend.justification = "left",
      legend.margin = margin(t = 0, r = 0, b = 4, l = 0),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 1, face = "bold"),
      strip.placement = "outside",
      plot.title = element_text(size = base_size + 3, face = "bold", 
                                hjust = 0, margin = margin(t = 0, r = 0, b = 4, l = 0)),
      plot.subtitle = element_text(size = base_size, color = "gray40", 
                                   hjust = 0, margin = margin(t = 0, r = 0, b = 6, l = 0)),
      plot.caption = element_text(size = base_size - 1, color = "gray50", 
                                  hjust = 0, margin = margin(t = 4, r = 0, b = 0, l = 0)),
      plot.margin = margin(t = 6, r = 10, b = 6, l = 6),
      complete = TRUE
    )
}

# ============================================
# 1.5 Earth-tone color scheme
# ============================================

#' Earth-tone palette
#' Inspired by natural colors: soil, vegetation, desert, rocks etc.
earth_palette <- function(n) {
  # Predefined earth-tone colors
  earth_colors <- c(
    "#8B5E3C",  # Brown (soil)
    "#C49A6C",  # Light brown (sand)
    "#4A7C59",  # Dark green (forest)
    "#7CAF6A",  # Medium green (vegetation)
    "#D4A574",  # Beige (desert)
    "#6B7B8D",  # Gray-blue (rock)
    "#A67B5B",  # Coffee (bark)
    "#5C6B4A",  # Olive green
    "#B89A7A",  # Khaki (dry grass)
    "#8B7D6B",  # Gray-brown (dead wood)
    "#3D5A3B",  # Ink green (deep forest)
    "#D4B896",  # Wheat
    "#7A6B5A",  # Gray-brown
    "#A0522D",  # Ochre
    "#6B8E6B"   # Sage green
  )
  
  if (n <= length(earth_colors)) {
    return(earth_colors[1:n])
  } else {
    # If more colors needed, use color space interpolation
    return(colorRampPalette(earth_colors)(n))
  }
}

#' Earth-tone color scheme (named colors)
earth_colors_named <- function() {
  c(
    # Two-group comparison
    "Low Zn" = "#4A7C59",      # Forest green
    "High Zn" = "#8B5E3C",     # Soil brown
    
    # Response variable categories
    "Biomass" = "#2D5A27",     # Dark green
    "Zn_concentration" = "#C49A6C",  # Sand
    "BCF" = "#6B7B8D",         # Gray-blue
    "Phosphorus" = "#A67B5B",  # Coffee
    
    # Direction
    "High > Low" = "#8B5E3C",  # Soil brown
    "Low > High" = "#4A7C59",  # Forest green
    "Not significant" = "#D4C5B0"  # Beige-gray
  )
}

#' Heatmap gradient (earth-tone)
#' From dark green to beige to dark brown
earth_gradient <- function() {
  c("#2D5A27", "#7CAF6A", "#E8D5B5", "#C49A6C", "#8B5E3C")
}

#' Significance stars
sig_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# ============================================
# 2. Statistical enhancement functions
# ============================================

#' Compute standard error from confidence interval (using t distribution)
get_se_from_ci <- function(lower, upper, ci_level = CI_LEVEL, df = NULL) {
  alpha <- 1 - ci_level
  if (!is.null(df) && df > 0 && is.finite(df)) {
    crit <- qt(1 - alpha/2, df = df)
  } else {
    crit <- qnorm(1 - alpha/2)
  }
  return((upper - lower) / (2 * crit))
}

#' Clean response variable names
clean_response <- function(resp) {
  recode(resp,
         "Biomass_shoot" = "Shoot biomass",
         "Biomass_root" = "Root biomass",
         "Biomass_fruit" = "Fruit biomass",
         "Biomass_leave" = "Leaf biomass",
         "Plant_height" = "Plant height",
         "Zn_shoot" = "Shoot Zn",
         "Zn_root" = "Root Zn",
         "Zn_plant" = "Plant Zn",
         "BCF_shoot" = "BCF (shoot)",
         "BCF_root" = "BCF (root)",
         "BCF_plant" = "BCF (plant)",
         "P_shoot" = "Shoot P",
         "P_root" = "Root P",
         .default = resp
  )
}

#' Two-level random-effects meta-regression comparison (core function)
compare_groups_bilevel <- function(low_data, high_data) {
  
  if (is.null(low_data) || is.null(high_data)) {
    return(NULL)
  }
  
  if (!is.null(low_data$effect) && !is.null(high_data$effect) &&
      !is.null(low_data$var) && !is.null(high_data$var)) {
    
    combined <- data.frame(
      effect = c(low_data$effect, high_data$effect),
      var = c(low_data$var, high_data$var),
      group = factor(c(rep("Low", length(low_data$effect)), 
                       rep("High", length(high_data$effect)))),
      study_id = c(low_data$study_id, high_data$study_id)
    )
    
    if (!is.null(low_data$effect_id) && !is.null(high_data$effect_id)) {
      combined$effect_id <- c(low_data$effect_id, high_data$effect_id)
    } else {
      combined$effect_id <- 1:nrow(combined)
    }
    
    combined <- combined[complete.cases(combined[, c("effect", "var", "study_id")]), ]
    
    if (nrow(combined) >= 4) {
      model <- tryCatch({
        rma.mv(
          yi = effect,
          V = var,
          mods = ~ group,
          random = list(~ 1 | study_id, ~ 1 | effect_id),
          data = combined,
          method = "REML",
          control = list(maxiter = 1000)
        )
      }, error = function(e) {
        tryCatch({
          rma(yi = effect, vi = var, mods = ~ group, data = combined, method = "REML")
        }, error = function(e2) NULL)
      })
      
      if (!is.null(model)) {
        is_multilevel <- inherits(model, "rma.mv")
        beta <- model$b[2]
        se <- model$se[2]
        p_val <- model$pval[2]
        
        if (is_multilevel && !is.null(model$dfs)) {
          df <- model$dfs[2]
        } else {
          df <- model$k - model$p - 1
        }
        
        crit <- qt(1 - (1 - CI_LEVEL)/2, df = max(1, df))
        ci_lower <- beta - crit * se
        ci_upper <- beta + crit * se
        
        if (is_multilevel) {
          tau2_total <- sum(model$sigma2)
        } else {
          tau2_total <- model$tau2
        }
        pi_lower <- beta - crit * sqrt(se^2 + tau2_total)
        pi_upper <- beta + crit * sqrt(se^2 + tau2_total)
        
        if (is_multilevel) {
          total_var <- sum(model$sigma2) + sum(model$tau2)
          i2 <- ifelse(total_var > 0, model$tau2[1] / total_var * 100, NA)
        } else {
          i2 <- model$I2
        }
        
        return(list(
          method = ifelse(is_multilevel, "Bilevel REML", "REML"),
          diff = as.numeric(beta),
          se = as.numeric(se),
          p_value = as.numeric(p_val),
          df = df,
          ci_lower = as.numeric(ci_lower),
          ci_upper = as.numeric(ci_upper),
          pi_lower = as.numeric(pi_lower),
          pi_upper = as.numeric(pi_upper),
          i2 = as.numeric(i2),
          k = model$k,
          n_low = length(unique(combined$study_id[combined$group == "Low"])),
          n_high = length(unique(combined$study_id[combined$group == "High"])),
          significant = p_val < 0.05,
          stars = sig_stars(p_val)
        ))
      }
    }
  }
  
  # Alternative method
  low_effect <- low_data$overall_pct
  high_effect <- high_data$overall_pct
  
  if (is.null(low_effect) || is.null(high_effect)) {
    return(NULL)
  }
  
  low_df <- if (!is.null(low_data$n_obs)) low_data$n_obs - 1 else 30
  high_df <- if (!is.null(high_data$n_obs)) high_data$n_obs - 1 else 30
  
  low_se <- get_se_from_ci(
    low_data$overall_pct_lower,
    low_data$overall_pct_upper,
    df = low_df
  )
  high_se <- get_se_from_ci(
    high_data$overall_pct_lower,
    high_data$overall_pct_upper,
    df = high_df
  )
  
  diff <- high_effect - low_effect
  diff_se <- sqrt(low_se^2 + high_se^2)
  
  df <- (low_se^2 + high_se^2)^2 / 
    (low_se^4 / low_df + high_se^4 / high_df)
  
  t_stat <- diff / diff_se
  p_value <- 2 * (1 - pt(abs(t_stat), df = df))
  
  crit <- qt(1 - (1 - CI_LEVEL)/2, df = df)
  ci_lower <- diff - crit * diff_se
  ci_upper <- diff + crit * diff_se
  
  pi_lower <- diff - qt(0.975, df = df) * diff_se * sqrt(1 + 1/df)
  pi_upper <- diff + qt(0.975, df = df) * diff_se * sqrt(1 + 1/df)
  
  return(list(
    method = "Welch t-test",
    diff = diff,
    se = diff_se,
    p_value = p_value,
    df = df,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    pi_lower = pi_lower,
    pi_upper = pi_upper,
    significant = p_value < 0.05,
    stars = sig_stars(p_value),
    t_stat = t_stat
  ))
}

# ============================================
# 3. Data loading and processing
# ============================================

safe_load_rds <- function(file_path, default = NULL) {
  if(file.exists(file_path)) {
    readRDS(file_path)
  } else {
    cat(sprintf("  ⚠️ File does not exist: %s\n", file_path))
    return(default)
  }
}

load_group_data <- function() {
  
  cat("\n========== Loading two-group data ==========\n")
  
  low_data <- safe_load_rds(file.path(PROCESSED_DIR, "overall_results_low.rds"))
  high_data <- safe_load_rds(file.path(PROCESSED_DIR, "overall_results_high.rds"))
  
  if (is.null(low_data) || is.null(high_data)) {
    cat("\n⚠️ Cannot load grouped data\n")
    return(NULL)
  }
  
  cat("\n✅ Data loading complete\n")
  cat(sprintf("   Low Zn:  %d response variables\n", length(low_data)))
  cat(sprintf("   High Zn: %d response variables\n", length(high_data)))
  
  return(list(low = low_data, high = high_data))
}

# ============================================
# 4. Prepare comparison data
# ============================================

prepare_comparison_data <- function(low_data, high_data) {
  
  cat("\n--- Preparing comparison data ---\n")
  
  common_resp <- intersect(names(low_data), names(high_data))
  cat(sprintf("  Common response variables: %d\n", length(common_resp)))
  
  comparison_list <- list()
  
  for(resp in common_resp) {
    
    low_res <- low_data[[resp]]
    high_res <- high_data[[resp]]
    
    if (is.null(low_res) || is.null(high_res)) {
      next
    }
    
    test_result <- compare_groups_bilevel(low_res, high_res)
    
    if (is.null(test_result)) {
      next
    }
    
    low_se <- get_se_from_ci(
      low_res$overall_pct_lower,
      low_res$overall_pct_upper,
      df = low_res$n_obs - 1
    )
    high_se <- get_se_from_ci(
      high_res$overall_pct_lower,
      high_res$overall_pct_upper,
      df = high_res$n_obs - 1
    )
    
    comparison_list[[resp]] <- data.frame(
      Response = resp,
      Category = ifelse(!is.null(low_res$category), low_res$category, "Unknown"),
      Low_Effect = low_res$overall_pct,
      Low_Lower = low_res$overall_pct_lower,
      Low_Upper = low_res$overall_pct_upper,
      Low_SE = low_se,
      Low_n = low_res$n_obs,
      Low_I2 = ifelse(!is.null(low_res$I2), low_res$I2, NA),
      High_Effect = high_res$overall_pct,
      High_Lower = high_res$overall_pct_lower,
      High_Upper = high_res$overall_pct_upper,
      High_SE = high_se,
      High_n = high_res$n_obs,
      High_I2 = ifelse(!is.null(high_res$I2), high_res$I2, NA),
      Diff = test_result$diff,
      Diff_SE = test_result$se,
      Diff_Lower = test_result$ci_lower,
      Diff_Upper = test_result$ci_upper,
      PI_Lower = test_result$pi_lower,
      PI_Upper = test_result$pi_upper,
      P_Value = test_result$p_value,
      Method = test_result$method,
      Significant = test_result$significant,
      Stars = test_result$stars,
      stringsAsFactors = FALSE
    )
  }
  
  comparison_data <- do.call(rbind, comparison_list)
  
  if (is.null(comparison_data) || nrow(comparison_data) == 0) {
    cat("⚠️ No data\n")
    return(NULL)
  }
  
  comparison_data$P_Adj <- p.adjust(comparison_data$P_Value, method = P_ADJ_METHOD)
  comparison_data$Significant_Adj <- comparison_data$P_Adj < 0.05
  comparison_data$Stars_Adj <- sig_stars(comparison_data$P_Adj)
  
  comparison_data$Category <- factor(comparison_data$Category, 
                                     levels = c("Biomass", "Zn_concentration", "BCF", "Phosphorus"))
  comparison_data$Response_clean <- clean_response(comparison_data$Response)
  
  comparison_data$Direction <- case_when(
    comparison_data$Diff > 0 & comparison_data$Significant_Adj ~ "High > Low",
    comparison_data$Diff < 0 & comparison_data$Significant_Adj ~ "Low > High",
    TRUE ~ "Not significant"
  )
  
  cat(sprintf("  ✅ Preparation complete: %d observations\n", nrow(comparison_data)))
  cat(sprintf("  📊 Significant differences (FDR): %d (%.1f%%)\n", 
              sum(comparison_data$Significant_Adj),
              sum(comparison_data$Significant_Adj) / nrow(comparison_data) * 100))
  
  return(comparison_data)
}

# ============================================
# 4.5 Fixed response variable order function
# ============================================

#' Get fixed response variable order
#' Fig1, Fig2, Fig4 all use this order
#' Note: this order is bottom-to-top (i.e., factor level order)
get_unified_response_order <- function(data) {
  
  # Define fixed order within each category (bottom to top)
  # This order applies to Fig1 y-axis (bottom to top) and Fig2/Fig4 y-axis (top to bottom)
  category_order <- list(
    "Biomass" = c(
      "Shoot biomass",      # Bottom
      "Fruit biomass",
      "Root biomass",
      "Plant height",
      "Leaf biomass"        # Top
    ),
    "Zn_concentration" = c(
      "Plant Zn",
      "Shoot Zn",
      "Root Zn"
    ),
    "BCF" = c(
      "BCF (plant)",
      "BCF (shoot)",
      "BCF (root)"
    ),
    "Phosphorus" = c(
      "Shoot P",
      "Root P"
    )
  )
  
  # Build full order: from Biomass to Phosphorus, within each category follow above order
  full_order <- c()
  for (cat in c("Biomass", "Zn_concentration", "BCF", "Phosphorus")) {
    if (cat %in% names(category_order)) {
      full_order <- c(full_order, category_order[[cat]])
    }
  }
  
  # Keep only response variables actually present in data
  existing_resp <- unique(data$Response_clean)
  full_order <- full_order[full_order %in% existing_resp]
  
  # If data has variables not in predefined order, append to end
  extra <- setdiff(existing_resp, full_order)
  if (length(extra) > 0) {
    cat(sprintf("  ⚠️ The following response variables are not in predefined order, appending to end: %s\n", 
                paste(extra, collapse = ", ")))
    full_order <- c(full_order, extra)
  }
  
  cat("\n--- Fixed response variable order (bottom to top) ---\n")
  for (i in seq_along(full_order)) {
    cat(sprintf("  %d. %s\n", i, full_order[i]))
  }
  
  return(full_order)
}

# ============================================
# 5. Fig1: Effect difference forest plot (earth-tone)
# ============================================

plot_forest <- function(data, response_order = NULL) {
  
  cat("\n--- Drawing effect difference forest plot (earth-tone) ---\n")
  
  # Use fixed order
  if (is.null(response_order)) {
    response_order <- get_unified_response_order(data)
  }
  
  # Set factor levels by fixed order (bottom to top)
  data$Response_clean <- factor(data$Response_clean, levels = response_order)
  
  # Ensure data sorted by factor levels
  data <- data[order(data$Response_clean), ]
  
  # Earth-tone colors
  sig_colors <- c(
    "High > Low" = "#8B5E3C",   # Soil brown
    "Low > High" = "#4A7C59",   # Forest green
    "Not significant" = "#D4C5B0"  # Beige-gray
  )
  
  category_positions <- aggregate(
    as.numeric(Response_clean) ~ Category, 
    data = data, 
    FUN = function(x) mean(x)
  )
  colnames(category_positions) <- c("Category", "mid")
  
  x_max <- max(abs(data$Diff_Upper), abs(data$Diff_Lower), na.rm = TRUE) * 1.3
  x_max <- max(x_max, 5)
  
  p <- ggplot(data, aes(x = Diff, y = Response_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    
    geom_errorbar(
      aes(xmin = PI_Lower, xmax = PI_Upper),
      width = 0.1, linewidth = 0.3, color = "gray70", alpha = 0.5,
      orientation = "y"
    ) +
    
    geom_errorbar(
      aes(xmin = Diff_Lower, xmax = Diff_Upper, color = Direction),
      width = 0.2, linewidth = 0.8, orientation = "y"
    ) +
    
    geom_point(
      aes(fill = Direction, size = Low_n + High_n),
      shape = 21, stroke = 0.6, color = "black", alpha = 0.9
    ) +
    
    geom_text(
      data = category_positions,
      aes(x = -x_max * 0.85, y = mid, label = Category),
      hjust = 0, size = 3.5, fontface = "bold", color = "black"
    ) +
    
    geom_text(
      data = data[data$Significant_Adj, ],
      aes(x = Diff_Upper + x_max * 0.03, label = Stars_Adj),
      size = 4, color = "#A0522D", fontface = "bold", hjust = 0
    ) +
    
    scale_color_manual(values = sig_colors, name = "Comparison", drop = FALSE) +
    scale_fill_manual(values = sig_colors, guide = "none", drop = FALSE) +
    
    scale_size_continuous(
      range = c(1.5, 5),
      name = "Sample size",
      breaks = c(10, 30, 50, 80),
      labels = c("10", "30", "50", "80")
    ) +
    
    scale_x_continuous(
      limits = c(-x_max, x_max),
      labels = function(x) paste0(x, "%")
    ) +
    
    labs(
      title = "AMF effect difference: High Zn vs Low Zn",
      subtitle = paste0("FDR-adjusted p-values | Prediction intervals (gray) | n = ", nrow(data), " comparisons"),
      x = "Effect difference (High Zn - Low Zn, %)",
      y = NULL,
      caption = "* p < 0.05, ** p < 0.01, *** p < 0.001 (FDR-adjusted)"
    ) +
    
    theme_journal(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      plot.margin = margin(t = 6, r = 20, b = 6, l = 6)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 3)),
      size = guide_legend(override.aes = list(shape = 21, fill = "gray70"))
    )
  
  return(p)
}

# ============================================
# 6. Fig2: Two-group effect size comparison (earth-tone)
# ============================================

plot_two_group <- function(data, response_order = NULL) {
  
  cat("\n--- Drawing two-group effect size comparison (earth-tone) ---\n")
  
  # Use fixed order
  if (is.null(response_order)) {
    response_order <- get_unified_response_order(data)
  }
  
  # Set factor levels by full order (top to bottom)
  data$Response_clean <- factor(data$Response_clean, levels = response_order)
  
  long_data <- data %>%
    dplyr::select(Response_clean, Category, Low_Effect, High_Effect,
                  Low_Lower, Low_Upper, High_Lower, High_Upper) %>%
    tidyr::pivot_longer(
      cols = c(Low_Effect, High_Effect),
      names_to = "Group",
      values_to = "Effect"
    ) %>%
    dplyr::mutate(
      Group = ifelse(Group == "Low_Effect", "Low Zn", "High Zn"),
      Group = factor(Group, levels = c("Low Zn", "High Zn")),
      Lower = ifelse(Group == "Low Zn", Low_Lower, High_Lower),
      Upper = ifelse(Group == "Low Zn", Low_Upper, High_Upper)
    )
  
  # Add category label positions for each response variable
  category_positions <- data %>%
    group_by(Category) %>%
    summarise(mid = mean(as.numeric(factor(Response_clean, levels = response_order))), .groups = "drop")
  
  # Earth-tone colors
  group_colors <- c("Low Zn" = "#4A7C59", "High Zn" = "#8B5E3C")
  
  # Category colors (for background or labels)
  category_colors <- c(
    "Biomass" = "#2D5A27",
    "Zn_concentration" = "#C49A6C",
    "BCF" = "#6B7B8D",
    "Phosphorus" = "#A67B5B"
  )
  
  p <- ggplot(long_data, aes(x = Effect, y = Response_clean, color = Group, shape = Group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_errorbar(
      aes(xmin = Lower, xmax = Upper),
      width = 0.2, linewidth = 0.4, alpha = 0.6, orientation = "y"
    ) +
    geom_point(size = 2.5, stroke = 0.3) +
    
    # Add category labels (left of y-axis)
    geom_text(
      data = category_positions,
      aes(x = -Inf, y = mid, label = Category),
      hjust = -0.1, size = 3.5, fontface = "bold", 
      color = "black", inherit.aes = FALSE
    ) +
    
    scale_color_manual(values = group_colors, name = "Zn level") +
    scale_shape_manual(values = c("Low Zn" = 16, "High Zn" = 17), name = "Zn level") +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    
    # Add category separator lines
    geom_hline(
      data = category_positions %>% filter(Category != "Phosphorus"),
      aes(yintercept = mid + 0.5),
      linetype = "dotted", color = "gray70", linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    
    labs(
      title = "AMF effects: Low Zn vs High Zn",
      subtitle = "Comparison across response variables",
      x = "Effect size (%)",
      y = NULL
    ) +
    theme_journal(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 8),
      panel.grid.major.x = element_line(color = "gray95", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      plot.margin = margin(t = 6, r = 10, b = 6, l = 30)
    )
  
  return(p)
}

# ============================================
# 7. Fig3: Scatter plot (earth-tone)
# ============================================

plot_scatter <- function(data) {
  
  cat("\n--- Drawing scatter plot (earth-tone) ---\n")
  
  # Earth-tone category colors
  category_colors <- c(
    "Biomass" = "#2D5A27",       # Dark green
    "Zn_concentration" = "#C49A6C",  # Sand
    "BCF" = "#6B7B8D",           # Gray-blue
    "Phosphorus" = "#A67B5B"     # Coffee
  )
  
  all_vals <- c(data$Low_Effect, data$High_Effect, 0)
  range_val <- diff(range(all_vals))
  range_val <- max(range_val, 10)
  lim <- c(min(all_vals) - range_val * 0.15, max(all_vals) + range_val * 0.15)
  
  p <- ggplot(data, aes(x = Low_Effect, y = High_Effect)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray80", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray80", linewidth = 0.3) +
    
    geom_point(
      aes(fill = Category, size = abs(Diff), shape = Significant_Adj),
      color = "black", stroke = 0.5, alpha = 0.85
    ) +
    
    geom_text_repel(
      aes(label = Response_clean, color = Category),
      size = 2.8, max.overlaps = 15, show.legend = FALSE,
      fontface = "plain", segment.size = 0.3
    ) +
    
    scale_fill_manual(values = category_colors, name = "Category") +
    scale_color_manual(values = category_colors, guide = "none") +
    
    scale_shape_manual(
      values = c(`TRUE` = 21, `FALSE` = 21),
      labels = c(`TRUE` = "Significant (FDR)", `FALSE` = "Not significant"),
      name = "Significance"
    ) +
    
    scale_size_continuous(
      range = c(1.5, 6),
      name = "|Difference| (%)",
      breaks = c(5, 15, 25, 35),
      labels = function(x) paste0(x, "%")
    ) +
    
    scale_x_continuous(limits = lim, labels = function(x) paste0(x, "%")) +
    scale_y_continuous(limits = lim, labels = function(x) paste0(x, "%")) +
    
    labs(
      title = "Low Zn vs High Zn AMF effects",
      subtitle = "Point size = absolute difference | Color = response category",
      x = "Effect in Low Zn (%)",
      y = "Effect in High Zn (%)",
      caption = "Points above diagonal: High Zn > Low Zn | Solid = FDR-significant"
    ) +
    
    theme_journal(base_size = 10) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      panel.grid.major = element_line(color = "gray95", linewidth = 0.3)
    )
  
  return(p)
}

# ============================================
# 7.5 Fig4: Difference heatmap (earth-tone) - corrected version
# ============================================

plot_heatmap_corrected <- function(data, response_order = NULL) {
  
  cat("\n--- Drawing difference heatmap (earth-tone) ---\n")
  
  if (is.null(response_order)) {
    response_order <- get_unified_response_order(data)
  }
  
  # Fig4: display top to bottom
  data$Response_clean <- factor(data$Response_clean, levels = response_order)
  
  # ==========================================
  # Key modification: change to High - Low
  # ==========================================
  data$Diff_High_vs_Low <- data$High_Effect - data$Low_Effect
  
  # Sort by factor levels
  data <- data[order(data$Response_clean), ]
  
  plot_data <- data %>%
    dplyr::select(Response_clean, Category, Diff_High_vs_Low, Significant_Adj) %>%
    dplyr::mutate(
      Label = ifelse(Significant_Adj, 
                     paste0(sprintf("%.1f", Diff_High_vs_Low), "*"),
                     sprintf("%.1f", Diff_High_vs_Low))
    )
  
  # Ensure Response_clean retains factor level order
  plot_data$Response_clean <- factor(plot_data$Response_clean, levels = response_order)
  plot_data$Category <- factor(plot_data$Category, 
                               levels = c("Biomass", "Zn_concentration", "BCF", "Phosphorus"))
  
  max_diff <- max(abs(plot_data$Diff_High_vs_Low), na.rm = TRUE)
  max_diff <- max(max_diff, 1)
  
  earth_gradient_colors <- c("#2D5A27", "#7CAF6A", "#E8D5B5", "#C49A6C", "#8B5E3C")
  
  p <- ggplot(plot_data, aes(x = "High - Low", y = Response_clean, fill = Diff_High_vs_Low)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = Label),
      size = 3,
      color = ifelse(abs(plot_data$Diff_High_vs_Low) > max_diff * 0.5, "white", "black")
    ) +
    scale_fill_gradientn(
      colors = earth_gradient_colors,
      limits = c(-max_diff, max_diff),
      name = "High Zn - Low Zn (%)"
    ) +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "AMF effect difference: High Zn vs Low Zn",
      subtitle = "* FDR-adjusted p < 0.05 | Brown: High Zn > Low Zn | Green: Low Zn > High Zn",
      x = NULL,
      y = NULL
    ) +
    theme_journal(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9, face = "bold"),
      strip.text = element_text(size = 9),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.key.height = unit(0.8, "cm"),
      panel.grid = element_blank()
    )
  
  return(p)
}

# ============================================
# 8. Main function
# ============================================

run_comparison <- function() {
  
  cat("\n========================================\n")
  cat("  Zn_Level two-group comparison analysis (earth-tone)\n")
  cat("  Compare AMF effect differences between Low Zn vs High Zn\n")
  cat("  Statistical method: two-level meta-regression + FDR correction + prediction intervals\n")
  cat("  Color scheme: Earth Tones\n")
  cat("========================================\n")
  
  group_data <- load_group_data()
  
  if (is.null(group_data)) {
    cat("\n⚠️ Cannot load grouped data\n")
    return(NULL)
  }
  
  comparison_data <- prepare_comparison_data(group_data$low, group_data$high)
  
  if (is.null(comparison_data)) {
    cat("\n⚠️ No data\n")
    return(NULL)
  }
  
  # ==========================================
  # Generate fixed response variable order
  # ==========================================
  response_order <- get_unified_response_order(comparison_data)
  
  # ==========================================
  # Path handling
  # ==========================================
  
  fig_dir <- file.path(getwd(), "results", "figures", "comparison")
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE)
    cat(sprintf("📁 Created directory: %s\n", fig_dir))
  }
  
  cat("\n--- Generating visualizations (earth-tone) ---\n")
  cat(sprintf("📁 Output directory: %s\n", fig_dir))
  
  # Fig1: Forest plot (using fixed order)
  p1 <- plot_forest(comparison_data, response_order)
  height1 <- max(4, nrow(comparison_data) * 0.35 + 1)
  
  tryCatch({
    ggsave(file.path(fig_dir, "01_forest.png"), p1, width = 8.5, height = height1, dpi = 400, bg = "white")
    ggsave(file.path(fig_dir, "01_forest.pdf"), p1, width = 8.5, height = height1, device = "pdf")
    cat("✅ Saved: 01_forest.png/pdf\n")
  }, error = function(e) {
    cat(sprintf("⚠️ Failed to save 01_forest: %s\n", e$message))
  })
  
  # Fig2: Two-group comparison (using fixed order, no reversal)
  p2 <- plot_two_group(comparison_data, response_order)
  height2 <- max(5, nrow(comparison_data) * 0.3 + 1.5)
  
  tryCatch({
    ggsave(file.path(fig_dir, "02_two_group.png"), p2, width = 8.5, height = height2, dpi = 400, bg = "white")
    ggsave(file.path(fig_dir, "02_two_group.pdf"), p2, width = 8.5, height = height2, device = "pdf")
    cat("✅ Saved: 02_two_group.png/pdf\n")
  }, error = function(e) {
    cat(sprintf("⚠️ Failed to save 02_two_group: %s\n", e$message))
  })
  
  # Fig3: Scatter plot
  p3 <- plot_scatter(comparison_data)
  
  tryCatch({
    ggsave(file.path(fig_dir, "03_scatter.png"), p3, width = 7, height = 6, dpi = 400, bg = "white")
    ggsave(file.path(fig_dir, "03_scatter.pdf"), p3, width = 7, height = 6, device = "pdf")
    cat("✅ Saved: 03_scatter.png/pdf\n")
  }, error = function(e) {
    cat(sprintf("⚠️ Failed to save 03_scatter: %s\n", e$message))
  })
  
  # Fig4: Heatmap (using fixed order, no reversal)
  p4 <- plot_heatmap_corrected(comparison_data, response_order)
  height4 <- max(5, nrow(comparison_data) * 0.25 + 1.5)
  
  tryCatch({
    ggsave(file.path(fig_dir, "04_heatmap.png"), p4, width = 7, height = height4, dpi = 400, bg = "white")
    ggsave(file.path(fig_dir, "04_heatmap.pdf"), p4, width = 7, height = height4, device = "pdf")
    cat("✅ Saved: 04_heatmap.png/pdf\n")
  }, error = function(e) {
    cat(sprintf("⚠️ Failed to save 04_heatmap: %s\n", e$message))
  })
  
  # Save data
  if(exists("SAVE_TABLES") && SAVE_TABLES) {
    tables_dir <- file.path(getwd(), "results", "tables")
    if(!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
    write.csv(comparison_data, file.path(tables_dir, "zn_level_comparison.csv"), row.names = FALSE)
    cat("\n✅ Saved data: zn_level_comparison.csv\n")
  }
  
  # Generate report
  cat("\n--- Generating report ---\n")
  
  report_lines <- c(
    "========================================",
    "Zn_Level two-group comparison analysis report (earth-tone)",
    paste("Generated at:", Sys.time()),
    "========================================",
    "",
    "1. Statistical methods:",
    paste("   - Primary method: two-level random-effects meta-regression (rma.mv)"),
    paste("   - Alternative method: Welch t-test"),
    paste("   - Multiple correction: FDR"),
    paste("   - Confidence interval level:", CI_LEVEL * 100, "%"),
    "",
    "2. Color scheme:",
    "   - Earth Tones",
    "   - Main colors: forest green (#4A7C59) and soil brown (#8B5E3C)",
    "   - Heatmap gradient: dark green -> beige -> dark brown",
    "",
    "3. Data summary:",
    paste("   - Number of response variables:", nrow(comparison_data)),
    paste("   - Category distribution:"),
    paste("     - Biomass:", sum(comparison_data$Category == "Biomass")),
    paste("     - Zn_concentration:", sum(comparison_data$Category == "Zn_concentration")),
    paste("     - BCF:", sum(comparison_data$Category == "BCF")),
    paste("     - Phosphorus:", sum(comparison_data$Category == "Phosphorus")),
    "",
    "4. Effect comparison:",
    paste("   - Low Zn group mean effect:", round(mean(comparison_data$Low_Effect, na.rm = TRUE), 1), "%"),
    paste("   - High Zn group mean effect:", round(mean(comparison_data$High_Effect, na.rm = TRUE), 1), "%"),
    paste("   - Mean difference (High - Low):", round(mean(comparison_data$Diff, na.rm = TRUE), 1), "%"),
    "",
    "5. Statistical significance (FDR-adjusted):",
    paste("   - Number of significant differences:", sum(comparison_data$Significant_Adj, na.rm = TRUE)),
    paste("   - High > Low:", sum(comparison_data$Diff > 0 & comparison_data$Significant_Adj, na.rm = TRUE)),
    paste("   - Low > High:", sum(comparison_data$Diff < 0 & comparison_data$Significant_Adj, na.rm = TRUE)),
    "",
    "6. Unified response variable order (top to bottom):",
    paste("   -", paste(response_order, collapse = " -> ")),
    "",
    "7. Visualizations:",
    "   - 01_forest.png/pdf (effect difference forest plot + prediction intervals)",
    "   - 02_two_group.png/pdf (two-group effect size comparison)",
    "   - 03_scatter.png/pdf (Low vs High scatter plot)",
    "   - 04_heatmap.png/pdf (High - Low difference heatmap)",
    "",
    "========================================"
  )
  
  report_path <- file.path(getwd(), "results", "zn_level_comparison_report.txt")
  report_dir <- dirname(report_path)
  if(!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)
  writeLines(report_lines, report_path)
  cat("✅ Report saved:", report_path, "\n")
  
  cat("\n========================================\n")
  cat("  ✅ Zn_Level two-group comparison analysis complete!\n")
  cat("========================================\n")
  cat("\n📁 Output directory:", fig_dir, "\n")
  cat("📄 Report:", report_path, "\n")
  
  invisible(comparison_data)
}

# ============================================
# 9. Run
# ============================================

if(exists("USE_ZN_LEVEL_GROUP") && USE_ZN_LEVEL_GROUP) {
  run_comparison()
} else {
  cat("\n⚠️ USE_ZN_LEVEL_GROUP = FALSE, skipping grouped comparison\n")
  cat("To run, set USE_ZN_LEVEL_GROUP <- TRUE in 00_config.R\n")
}