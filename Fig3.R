# ============================================
# Final integrated version | Three-level meta rma.mv + plant part comparison (quality-improved version)
# Improvements: de-duplicated statistical labels, unified theme, raw data visualization, enhanced error handling
# ============================================

# ============================================
# 1. Setup and package loading
# ============================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(metafor)
  library(patchwork)
  library(viridis)
  library(gridExtra)
  library(grid)
  library(stringr)
})

# Explicitly handle common naming conflicts
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

set.seed(20240729)

FIG_OUT_PATH <- "output_figures"
dir.create(FIG_OUT_PATH, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))

# Unified color palette (removed inconsistent yellow, All unified as gray)
GROUP_COLORS <- c("Low Zn" = "#440154", "High Zn" = "#21908C", "All" = "#757575")
GROUP_SHAPES <- c("Low Zn" = 16, "High Zn" = 17, "All" = 15)

# Unified theme function
theme_meta <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10, color = "gray20"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(t = 6),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.9, "lines"),
      strip.background = element_rect(fill = "#F5F5F5", color = NA),
      strip.text = element_text(size = 10, face = "bold"),
      plot.caption = element_text(hjust = 0, size = 8, color = "gray50", face = "italic")
    )
}

# ============================================
# 2. Basic helper functions
# ============================================
sig_star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***", 
                ifelse(p < 0.01, "**", 
                       ifelse(p < 0.05, "*", "ns"))))
}

format_p <- function(p_val) {
  if(is.na(p_val)) return("NA")
  if(p_val < 0.001) return("<0.001***")
  if(p_val < 0.01) return(paste0(round(p_val, 3), "**"))
  if(p_val < 0.05) return(paste0(round(p_val, 3), "*"))
  return(paste0(round(p_val, 3), "ns"))
}

fdr_correct <- function(p_vector){
  if(length(p_vector) == 0) return(numeric(0))
  na_idx <- is.na(p_vector)
  if(all(na_idx)) return(rep(NA, length(p_vector)))
  p_clean <- p_vector[!na_idx]
  p_adj <- p.adjust(p_clean, method = "BH")
  result <- rep(NA, length(p_vector))
  result[!na_idx] <- p_adj
  return(result)
}

safe_save <- function(plot_obj, filepath, w, h) {
  tryCatch({
    if (grepl("\\.png$", filepath, ignore.case = TRUE)) {
      ggsave(filepath, plot = plot_obj, width = w, height = h, dpi = 300, bg = "white")
    } else {
      ggsave(filepath, plot = plot_obj, width = w, height = h)
    }
    cat("✅ Saved:", filepath, "\n")
  }, error = function(e) {
    message("Save failed: ", filepath, " - ", e$message)
  })
}

# ============================================
# 3. Statistical functions
# ============================================
extract_coef <- function(mod, model_name, apply_fdr = TRUE){
  if(is.null(mod) || !inherits(mod, "rma.mv")) {
    return(data.frame(Model = model_name, Term = NA, Estimate = NA, SE = NA, Z = NA, 
                      pval = NA, p_fdr = NA, Significance = NA))
  }
  
  Term <- rownames(mod$b)
  Estimate <- as.vector(mod$b)
  SE <- mod$se
  Z <- mod$zval
  pval <- mod$pval
  
  est_df <- data.frame(
    Model = model_name,
    Term = Term,
    Estimate = Estimate,
    SE = SE,
    Z = Z,
    pval = pval,
    stringsAsFactors = FALSE
  )
  
  if(apply_fdr) {
    est_df$p_fdr <- fdr_correct(est_df$pval)
    est_df$Significance <- sig_star(est_df$p_fdr)
  } else {
    est_df$p_fdr <- est_df$pval
    est_df$Significance <- sig_star(est_df$pval)
  }
  
  return(est_df)
}

report_heterogeneity_fixed <- function(mod, mod_name = "Model") {
  if(is.null(mod) || !inherits(mod, "rma.mv")) {
    return(data.frame(
      Model = mod_name,
      tau2_study = NA, tau2_obs = NA,
      I2_study = NA, I2_obs = NA, I2_total = NA,
      mean_sampling_var = NA, total_variance = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  tryCatch({
    tau2_study <- mod$sigma2[1]
    tau2_obs   <- mod$sigma2[2]
    mean_vi    <- mean(mod$vi, na.rm = TRUE)
    total_var  <- tau2_study + tau2_obs + mean_vi
    
    if(total_var > 1e-9 && is.finite(total_var)){
      I2_study <- 100 * tau2_study / total_var
      I2_obs   <- 100 * tau2_obs / total_var
      I2_total <- 100 * (tau2_study + tau2_obs) / total_var
    } else {
      I2_study <- I2_obs <- I2_total <- 0
    }
    
    data.frame(
      Model = mod_name,
      tau2_study = round(tau2_study, 4),
      tau2_obs = round(tau2_obs, 4),
      I2_study = round(I2_study, 1),
      I2_obs = round(I2_obs, 1),
      I2_total = round(I2_total, 1),
      mean_sampling_var = round(mean_vi, 4),
      total_variance = round(total_var, 4),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("【Heterogeneity error - ", mod_name, "】", e$message)
    report_heterogeneity_fixed(NULL, mod_name)
  })
}

predict_meta_line_with_p <- function(dat, x_var, y_var, group_name = "") {
  if(nrow(dat) < 5) {
    message("⚠️ ", group_name, " insufficient data (n=", nrow(dat), "), skipping ", x_var)
    return(list(line = NULL, p_value = NA, significant = FALSE, slope = NA, n = nrow(dat)))
  }
  
  var_col <- case_when(
    y_var == "lnRR_biomass" ~ "variance_biomass",
    y_var == "lnRR_Zn" ~ "variance_Zn",
    y_var == "lnRR_P" ~ "variance_P",
    TRUE ~ "variance_biomass"
  )
  
  if(!var_col %in% colnames(dat)) {
    message("⚠️ ", group_name, " variance column ", var_col, " missing, skipping ", x_var)
    return(list(line = NULL, p_value = NA, significant = FALSE, slope = NA, n = nrow(dat)))
  }
  
  tryCatch({
    x_range <- range(dat[[x_var]], na.rm = TRUE)
    if(diff(x_range) < 0.001) {
      message("⚠️ ", group_name, " predictor ", x_var, " has near-zero variance, skipping")
      return(list(line = NULL, p_value = NA, significant = FALSE, slope = NA, n = nrow(dat)))
    }
    
    mod <- rma.mv(
      yi = dat[[y_var]],
      V = dat[[var_col]],
      mods = as.formula(paste0("~", x_var)),
      random = ~ 1 | study_id / obs_id, 
      data = dat, 
      method = "REML",
      control = list(optimizer = "bobyqa", maxiter = 10000, stepadj = 0.5)
    )
    
    if(!is.null(mod$converged) && !mod$converged) {
      message("⚠️ ", group_name, " model did not converge")
      return(list(line = NULL, p_value = NA, significant = FALSE, slope = NA, n = nrow(dat)))
    }
    
    p_val <- mod$pval[2]
    if(is.na(p_val) || length(p_val) == 0) p_val <- 1
    
    newx <- seq(x_range[1], x_range[2], length.out = 80)
    pred <- predict(mod, newmods = cbind(newx))
    
    line_df <- tibble(
      !!x_var := newx, 
      pred_lnRR = pred$pred, 
      ci_lb = pred$ci.lb, 
      ci_ub = pred$ci.ub,
      Group = group_name
    )
    
    return(list(
      line = line_df, 
      p_value = p_val, 
      significant = p_val < 0.05,
      slope = mod$b[2],
      n = nrow(dat)
    ))
    
  }, error = function(e) {
    message("⚠️ ", group_name, " fitting failed: ", e$message)
    return(list(line = NULL, p_value = NA, significant = FALSE, slope = NA, n = nrow(dat)))
  })
}

calc_mean_effect <- function(dat, resp_name, groupname){
  if(nrow(dat) < 5) {
    message("Warning: ", resp_name, " ", groupname, " only has ", nrow(dat), " observations")
    return(NULL)
  }
  tryCatch({
    mod <- rma.mv(yi = lnRR, V = variance, random = ~ 1 | study_id / obs_id, 
                  data = dat, method = "REML",
                  control = list(optimizer = "bobyqa", maxiter = 10000))
    pred <- predict(mod)
    tibble(Group = groupname, 
           Response = resp_name, 
           mean_lnRR = pred$pred, 
           ci_lb = pred$ci.lb, 
           ci_ub = pred$ci.ub, 
           p = mod$pval,
           n = nrow(dat))
  }, error = function(e) {
    message("Mean effect fitting failed: ", resp_name, " ", groupname, " - ", e$message)
    return(NULL)
  })
}

diagnose_subgroup <- function(mod, data, name) {
  cat("\n========== ", name, " Diagnostics ==========\n")
  
  converged <- if(!is.null(mod$converged)) mod$converged else NA
  cat("【Convergence】:", ifelse(isTRUE(converged), "✅ Converged", "❌ Not converged"), "\n")
  
  tryCatch({
    inf <- influence(mod)
    cook_cutoff <- 4 / nrow(data)
    high_cook <- which(inf$cook.d > cook_cutoff)
    cat("【High influence points】: ", length(high_cook), ", Cook's D threshold = ", round(cook_cutoff, 4), "\n")
  }, error = function(e) {
    cat("【High influence points】: Calculation failed (rma.mv influence limitation)\n")
  })
  
  tryCatch({
    cor_matrix <- cor(data[, c("lnRR_Zn", "lnRR_P")], use = "complete.obs")
    cat("【lnRR_Zn ↔ lnRR_P correlation r =", round(cor_matrix[1,2], 3), "]\n")
    if(abs(cor_matrix[1,2]) > 0.7){
      message("⚠️ Warning: High collinearity between predictors!")
    }
  }, error = function(e) {
    cat("【Variable correlation】Calculation failed\n")
  })
}

# ============================================
# 4. Plotting functions (improved version)
# ============================================
create_table_plot <- function(table_df, title, font_size = 4){
  tbl <- tableGrob(table_df, rows = NULL, 
                   theme = ttheme_default(base_size = font_size, 
                                          padding = unit(c(4, 4), "mm")))
  p <- ggplot() + 
    annotation_custom(tbl) +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
  return(p)
}

create_scatter_with_stats <- function(data, xcol, ycol, xlab, ylab, title) {
  low_dat <- data %>% filter(Zn_Level == "Low Zn")
  high_dat <- data %>% filter(Zn_Level == "High Zn")
  all_dat <- data
  
  res_all <- predict_meta_line_with_p(all_dat, xcol, ycol, "All")
  res_low <- predict_meta_line_with_p(low_dat, xcol, ycol, "Low Zn")
  res_high <- predict_meta_line_with_p(high_dat, xcol, ycol, "High Zn")
  
  n_low <- nrow(low_dat)
  n_high <- nrow(high_dat)
  n_all <- nrow(all_dat)
  
  status_all <- ifelse(is.null(res_all$line), "fail", ifelse(res_all$significant, "sig", "ns"))
  status_low <- ifelse(is.null(res_low$line), "fail", ifelse(res_low$significant, "sig", "ns"))
  status_high <- ifelse(is.null(res_high$line), "fail", ifelse(res_high$significant, "sig", "ns"))
  
  slope_all <- ifelse(is.null(res_all$slope) || is.na(res_all$slope), NA, res_all$slope)
  slope_low <- ifelse(is.null(res_low$slope) || is.na(res_low$slope), NA, res_low$slope)
  slope_high <- ifelse(is.null(res_high$slope) || is.na(res_high$slope), NA, res_high$slope)
  
  p_slope_all <- ifelse(is.null(res_all$p_value) || is.na(res_all$p_value), "NA", format_p(res_all$p_value))
  p_slope_low <- ifelse(is.null(res_low$p_value) || is.na(res_low$p_value), "NA", format_p(res_low$p_value))
  p_slope_high <- ifelse(is.null(res_high$p_value) || is.na(res_high$p_value), "NA", format_p(res_high$p_value))
  
  annotation_text <- paste0(
    "All    : β = ", round(slope_all, 3), ", p = ", p_slope_all, ", n = ", n_all, " [", status_all, "]\n",
    "Low Zn : β = ", round(slope_low, 3), ", p = ", p_slope_low, ", n = ", n_low, " [", status_low, "]\n",
    "High Zn: β = ", round(slope_high, 3), ", p = ", p_slope_high, ", n = ", n_high, " [", status_high, "]"
  )
  
  var_col <- case_when(
    ycol == "lnRR_biomass" ~ "variance_biomass",
    ycol == "lnRR_Zn" ~ "variance_Zn",
    ycol == "lnRR_P" ~ "variance_P",
    TRUE ~ "variance_biomass"
  )
  
  # Safely compute weights
  data_with_weight <- data %>%
    mutate(
      weight = 1 / .data[[var_col]],
      weight = ifelse(is.finite(weight) & weight > 0, weight, NA)
    ) %>%
    filter(!is.na(weight)) %>%
    mutate(
      weight_norm = 1 + 4 * (weight - min(weight)) / (max(weight) - min(weight) + 1e-10)
    )
  
  # Use quantiles to generate legend breaks (more robust)
  weight_quants <- quantile(data_with_weight$weight, probs = c(0.1, 0.5, 0.9))
  w_min <- min(data_with_weight$weight)
  w_max <- max(data_with_weight$weight)
  weight_breaks_norm <- 1 + 4 * (weight_quants - w_min) / (w_max - w_min + 1e-10)
  weight_breaks_norm <- pmax(1, pmin(5, weight_breaks_norm))
  weight_labels <- format(round(weight_quants, 0), big.mark = ",")
  
  x_range <- range(data[[xcol]], na.rm = TRUE)
  y_range <- range(data[[ycol]], na.rm = TRUE)
  x_pos <- x_range[1] + 0.02 * diff(x_range)
  y_pos <- y_range[2] - 0.02 * diff(y_range)
  
  line_types <- c("All" = "dashed", "Low Zn" = "solid", "High Zn" = "solid")
  
  p <- ggplot(data_with_weight, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.4) +
    geom_rug(aes(color = Zn_Level), sides = "bl", alpha = 0.25, length = unit(0.018, "npc")) +
    geom_point(aes(color = Zn_Level, shape = Zn_Level, size = weight_norm),
               alpha = 0.5, stroke = 0.8) +
    scale_color_manual(values = GROUP_COLORS, name = "Zn Level") +
    scale_shape_manual(values = GROUP_SHAPES, name = "Zn Level") +
    scale_size_continuous(
      range = c(1.5, 6),
      name = "Weight (1/variance)",
      breaks = weight_breaks_norm,
      labels = weight_labels
    )
  
  if(!is.null(res_all$line)){
    p <- p + 
      geom_ribbon(data = res_all$line, aes(x = .data[[xcol]], ymin = ci_lb, ymax = ci_ub),
                  fill = "gray50", alpha = 0.1, inherit.aes = FALSE) +
      geom_line(data = res_all$line, aes(x = .data[[xcol]], y = pred_lnRR, linetype = "All"),
                color = "gray30", linewidth = 1, alpha = 0.8, inherit.aes = FALSE)
  }
  if(!is.null(res_low$line)){
    p <- p + 
      geom_ribbon(data = res_low$line, aes(x = .data[[xcol]], ymin = ci_lb, ymax = ci_ub),
                  fill = GROUP_COLORS["Low Zn"], alpha = 0.12, inherit.aes = FALSE) +
      geom_line(data = res_low$line, aes(x = .data[[xcol]], y = pred_lnRR, linetype = "Low Zn"),
                color = GROUP_COLORS["Low Zn"], linewidth = 1, alpha = 0.8, inherit.aes = FALSE)
  }
  if(!is.null(res_high$line)){
    p <- p + 
      geom_ribbon(data = res_high$line, aes(x = .data[[xcol]], ymin = ci_lb, ymax = ci_ub),
                  fill = GROUP_COLORS["High Zn"], alpha = 0.12, inherit.aes = FALSE) +
      geom_line(data = res_high$line, aes(x = .data[[xcol]], y = pred_lnRR, linetype = "High Zn"),
                color = GROUP_COLORS["High Zn"], linewidth = 1, alpha = 0.8, inherit.aes = FALSE)
  }
  
  p + 
    scale_linetype_manual(values = line_types, name = "Fitted line") +
    annotate("text", x = x_pos, y = y_pos,
             label = annotation_text, hjust = 0, vjust = 1, 
             size = 3, fontface = "bold", color = "gray20") +
    labs(x = xlab, y = ylab, title = title,
         caption = "Point size ∝ weight (1/variance) | Solid = p<0.05, Dashed = ns | β = meta-regression slope") +
    theme_meta() +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 3)),
      shape = guide_legend(order = 1, override.aes = list(size = 3)),
      linetype = guide_legend(order = 2),
      size = guide_legend(order = 3)
    )
}

# ============================================
# Figure 4 related functions (optimized version)
# ============================================
get_part_within_p <- function(data, response_var, variance_var, part) {
  sub_data <- data %>% filter(plant_part == part)
  if(nrow(sub_data) < 3) return(NA)
  tryCatch({
    mod <- rma.mv(
      yi = sub_data[[response_var]],
      V = sub_data[[variance_var]],
      random = ~ 1 | study_id / obs_id,
      data = sub_data,
      method = "REML",
      control = list(optimizer = "bobyqa", maxiter = 10000)
    )
    return(mod$pval)
  }, error = function(e) return(NA))
}

prepare_part_data <- function(all_data, response_var) {
  
  variance_var <- paste0("variance_", gsub("lnRR_", "", response_var))
  
  part_data <- all_data %>%
    filter(!is.na(.data[[response_var]]),
           !is.na(.data[[variance_var]]),
           !is.na(plant_part),          # Ensure no NA in plant part
           !is.na(Zn_Level),            # Ensure no NA in Zn level
           plant_part %in% c("root", "shoot"))
  
  if(nrow(part_data) < 10) return(NULL)
  
  # Create All group (merged only from Low Zn and High Zn, no NA)
  all_group <- part_data %>%
    mutate(Zn_Level = "All")
  
  full_data <- bind_rows(part_data, all_group) %>%
    mutate(
      Zn_Level = factor(Zn_Level, levels = c("Low Zn", "High Zn", "All")),
      plant_part = factor(plant_part, levels = c("root", "shoot"))
    )
  
  part_summary <- full_data %>%
    group_by(Zn_Level, plant_part) %>%
    summarise(
      n = n(),
      mean = mean(.data[[response_var]], na.rm = TRUE),
      se = sd(.data[[response_var]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Compute p-values for comparison against 0
  p_values <- full_data %>%
    group_by(Zn_Level, plant_part) %>%
    summarise(
      p = get_part_within_p(
        data.frame(
          lnRR = .data[[response_var]],
          variance = .data[[variance_var]],
          study_id = study_id,
          obs_id = obs_id,
          plant_part = plant_part
        ),
        "lnRR", "variance", plant_part[1]
      ),
      .groups = "drop"
    )
  
  part_summary <- left_join(part_summary, p_values, by = c("Zn_Level", "plant_part"))
  part_summary$p_fdr <- fdr_correct(part_summary$p)
  part_summary$Sig <- sig_star(part_summary$p_fdr)
  
  return(part_summary)
}

# 【Fixed】Merge ANOVA and pairwise comparison into a single Part effect
calculate_part_stats <- function(all_data, response_var) {
  variance_var <- paste0("variance_", gsub("lnRR_", "", response_var))
  
  part_data <- all_data %>%
    filter(!is.na(.data[[response_var]]),
           !is.na(.data[[variance_var]]),
           !is.na(plant_part),
           plant_part %in% c("root", "shoot"))
  
  if(nrow(part_data) < 10) return(NULL)
  
  all_group <- part_data %>% mutate(Zn_Level = "All")
  full_data <- bind_rows(part_data, all_group) %>%
    mutate(Zn_Level = factor(Zn_Level, levels = c("Low Zn", "High Zn", "All")))
  
  stats_results <- data.frame()
  for(g in levels(full_data$Zn_Level)) {
    sub_data <- full_data %>% filter(Zn_Level == g)
    
    part_effect <- tryCatch({
      mod <- rma.mv(
        yi = sub_data[[response_var]],
        V = sub_data[[variance_var]],
        mods = ~ plant_part,
        random = ~ 1 | study_id / obs_id,
        data = sub_data,
        method = "REML",
        control = list(optimizer = "bobyqa", maxiter = 10000)
      )
      p_val <- if(length(mod$pval) >= 2) mod$pval[2] else NA
      estimate <- if(length(mod$b) >= 2) mod$b[2] else NA
      se <- if(length(mod$se) >= 2) mod$se[2] else NA
      list(p = p_val, estimate = estimate, se = se)
    }, error = function(e) list(p = NA, estimate = NA, se = NA))
    
    stats_results <- bind_rows(stats_results, data.frame(
      Zn_Level = g, 
      part_p = part_effect$p,
      part_estimate = part_effect$estimate,
      part_se = part_effect$se,
      stringsAsFactors = FALSE
    ))
  }
  stats_results <- stats_results %>%
    mutate(part_p_fdr = fdr_correct(part_p))
  
  stats_results %>%
    mutate(
      Zn_Level = factor(Zn_Level, levels = c("Low Zn", "High Zn", "All")),
      part_label = case_when(
        is.na(part_p_fdr) ~ "ns",
        part_p_fdr < 0.001 ~ "***",
        part_p_fdr < 0.01 ~ "**",
        part_p_fdr < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      part_text = ifelse(is.na(part_p_fdr), 
                         "Part effect: ns",
                         sprintf("Part effect: %s (p = %.3f)", part_label, part_p_fdr))
    )
}

# ============================================
# Figure 4A: Biomass plant part comparison (no scatter, no NA group)
# ============================================
create_part_plot_biomass <- function(all_data) {
  
  # Explicitly filter NA to ensure they don't participate in analysis or plotting
  part_summary <- prepare_part_data(
    all_data %>% filter(!is.na(Zn_Level), !is.na(plant_part)), 
    "lnRR_biomass"
  )
  
  if(is.null(part_summary)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") + theme_void())
  }
  
  # Ensure summary table has no NA group again
  part_summary <- part_summary %>% filter(!is.na(Zn_Level), !is.na(plant_part))
  
  stat_labels <- calculate_part_stats(
    all_data %>% filter(!is.na(Zn_Level), !is.na(plant_part)), 
    "lnRR_biomass"
  )
  
  ymax <- max(abs(part_summary$mean) + part_summary$se, na.rm = TRUE) * 1.6
  if(is.na(ymax) || ymax < 0.1) ymax <- 1
  
  part_labels <- c("root" = "Root", "shoot" = "Shoot")
  
  p <- ggplot(part_summary, aes(x = plant_part, y = mean, fill = Zn_Level)) +
    geom_col(position = position_dodge(0.8), width = 0.7, 
             color = "gray30", linewidth = 0.3, alpha = 0.9) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(0.8), 
                  width = 0.15, linewidth = 0.7) +
    geom_text(aes(label = Sig, 
                  y = ifelse(mean >= 0, mean + se + 0.04 * ymax, mean - se - 0.04 * ymax)),
              position = position_dodge(0.8),
              size = 5, fontface = "bold", 
              vjust = ifelse(part_summary$mean >= 0, 0, 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    scale_fill_manual(values = GROUP_COLORS, name = "Zn Level") +
    scale_x_discrete(labels = part_labels) +
    labs(title = "A. Biomass Effect by Plant Part",
         x = "Plant Part", 
         y = expression("Mean lnRR"[Biomass] ~ "± SE")) +
    coord_cartesian(ylim = c(-ymax * 0.15, ymax)) +
    theme_meta() +
    theme(panel.grid.major.x = element_blank(),
          plot.margin = margin(t = 40, r = 10, b = 10, l = 10))
  
  # Add plant part effect labels
  if(!is.null(stat_labels) && nrow(stat_labels) > 0) {
    stat_labels <- stat_labels %>% 
      filter(!is.na(Zn_Level), !is.na(part_p))
    
    if(nrow(stat_labels) > 0) {
      p <- p +
        geom_text(data = stat_labels,
                  aes(x = 1.5, y = ymax * 0.94, label = part_text),
                  hjust = 0.5, vjust = 1, 
                  size = 3.5, fontface = "bold", color = "#C0392B",
                  inherit.aes = FALSE)
    }
  }
  
  return(p)
}

# ============================================
# Figure 4B: Zn and P plant part comparison (no scatter, no NA group)
# ============================================
create_part_plot_znp <- function(all_data) {
  
  # Filter NA before passing
  clean_data <- all_data %>% filter(!is.na(Zn_Level), !is.na(plant_part))
  
  zn_summary <- prepare_part_data(clean_data, "lnRR_Zn")
  if(!is.null(zn_summary)) zn_summary$Variable <- "Zn"
  
  p_summary <- prepare_part_data(clean_data, "lnRR_P")
  if(!is.null(p_summary)) p_summary$Variable <- "P"
  
  if(is.null(zn_summary) && is.null(p_summary)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") + theme_void())
  }
  
  part_summary <- bind_rows(zn_summary, p_summary)
  
  # Ensure no NA rows
  part_summary <- part_summary %>% 
    filter(!is.na(Zn_Level), !is.na(plant_part), !is.na(Variable))
  
  if(nrow(part_summary) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") + theme_void())
  }
  
  ymax <- max(abs(part_summary$mean) + part_summary$se, na.rm = TRUE) * 1.5
  if(is.na(ymax) || ymax < 0.1) ymax <- 1
  
  var_colors <- c("Zn" = "#4DBBD5", "P" = "#E64B35")
  part_labels <- c("root" = "Root", "shoot" = "Shoot")
  
  p <- ggplot(part_summary, aes(x = plant_part, y = mean, fill = Variable)) +
    geom_col(position = position_dodge(0.8), width = 0.7,
             color = "gray30", linewidth = 0.3, alpha = 0.9) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(0.8),
                  width = 0.15, linewidth = 0.7) +
    geom_text(aes(label = ifelse(is.na(Sig), "", Sig),
                  y = ifelse(mean >= 0, mean + se + 0.04 * ymax, mean - se - 0.04 * ymax)),
              position = position_dodge(0.8),
              size = 4.5, fontface = "bold",
              vjust = ifelse(part_summary$mean >= 0, 0, 1),
              na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    scale_fill_manual(values = var_colors, name = "Variable") +
    scale_x_discrete(labels = part_labels) +
    facet_wrap(~ Zn_Level, nrow = 1, 
               labeller = labeller(Zn_Level = function(x) paste("Zn Level:", x))) +
    labs(title = "B. Zn and P Effects by Plant Part",
         x = "Plant Part", 
         y = expression("Mean lnRR ± SE")) +
    coord_cartesian(ylim = c(-ymax * 0.15, ymax)) +
    theme_meta() +
    theme(panel.grid.major.x = element_blank(),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
  
  return(p)
}
create_figure4_optimized <- function(all_data) {
  cat("\n========== Generating optimized Figure 4 ==========\n")
  
  if(!"plant_part" %in% colnames(all_data)) {
    cat("⚠️ Missing 'plant_part' column\n")
    return(NULL)
  }
  
  part_types <- unique(na.omit(all_data$plant_part))
  if(length(part_types) < 2) {
    cat("⚠️ Less than 2 plant part types detected\n")
    return(NULL)
  }
  
  cat("Detected parts:", paste(part_types, collapse = ", "), "\n")
  cat("Analyzing root and shoot only\n")
  
  # ---- [NEW] Export Fig 4 sample sizes (n) for each category to CSV ----
  clean_data_for_n <- all_data %>% filter(!is.na(Zn_Level), !is.na(plant_part))
  fig4_n_list <- list()
  
  n_biomass <- prepare_part_data(clean_data_for_n, "lnRR_biomass")
  if(!is.null(n_biomass)) {
    fig4_n_list[["Fig4A_biomass"]] <- n_biomass %>%
      select(Zn_Level, plant_part, n, mean, se, p, p_fdr, Sig) %>%
      mutate(Panel = "Fig4A_biomass", Response = "lnRR_biomass")
  }
  
  n_zn <- prepare_part_data(clean_data_for_n, "lnRR_Zn")
  if(!is.null(n_zn)) {
    fig4_n_list[["Fig4B_Zn"]] <- n_zn %>%
      select(Zn_Level, plant_part, n, mean, se, p, p_fdr, Sig) %>%
      mutate(Panel = "Fig4B_Zn", Response = "lnRR_Zn")
  }
  
  n_p <- prepare_part_data(clean_data_for_n, "lnRR_P")
  if(!is.null(n_p)) {
    fig4_n_list[["Fig4B_P"]] <- n_p %>%
      select(Zn_Level, plant_part, n, mean, se, p, p_fdr, Sig) %>%
      mutate(Panel = "Fig4B_P", Response = "lnRR_P")
  }
  
  if(length(fig4_n_list) > 0) {
    fig4_n_all <- bind_rows(fig4_n_list) %>%
      relocate(Panel, Response, Zn_Level, plant_part, n)
    write.csv(fig4_n_all,
              file.path(FIG_OUT_PATH, "B_Figure4_n_and_summary.csv"),
              row.names = FALSE)
    cat("\n✅ Fig 4 n summary saved: B_Figure4_n_and_summary.csv\n")
    print(fig4_n_all)
  }
  # ---- [NEW END] ----
  
  p4a <- create_part_plot_biomass(all_data)
  p4b <- create_part_plot_znp(all_data)
  
  Figure4 <- (p4a | p4b) +
    plot_layout(guides = "collect", widths = c(1, 1.3)) &
    theme(legend.position = "bottom")
  
  Figure4 <- Figure4 +
    plot_annotation(
      title = "Figure 4. Plant Part Comparison Across Zn Groups",
      subtitle = "Root vs Shoot differences (***p<0.001, **p<0.01, *p<0.05, ns = not significant)",
      theme = theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40")
      )
    )
  
  return(Figure4)
}

# ============================================
# 5. Data import and preprocessing
# ============================================
if(!file.exists("Three_InRR2.csv")) {
  stop("❌ Data file 'Three_InRR2.csv' not found in working directory: ", getwd())
}

data_clean <- read.csv("Three_InRR2.csv")

all_data <- data_clean %>%
  mutate(
    Zn_Level = factor(Zn_Level, levels = c("low", "high"), labels = c("Low Zn", "High Zn")),
    obs_id = row_number()
  )

has_plant_part <- "plant_part" %in% colnames(all_data)

if(has_plant_part) {
  all_data <- all_data %>%
    mutate(plant_part = tolower(trimws(plant_part))) %>%
    filter(!is.na(plant_part))
}

reg_data <- all_data %>%
  filter(!is.na(lnRR_biomass) & !is.na(variance_biomass),
         !is.na(lnRR_Zn) & !is.na(variance_Zn),
         !is.na(lnRR_P) & !is.na(variance_P))

cat("\n========== Dataset Info ==========\n")
cat("Total observations n =", nrow(reg_data), "\n")
cat("Unique studies k =", length(unique(reg_data$study_id)), "\n")
print(table(reg_data$Zn_Level))

low_dat <- reg_data %>% filter(Zn_Level == "Low Zn")
high_dat <- reg_data %>% filter(Zn_Level == "High Zn")
all_dat <- reg_data

cat("\nSubgroup sample sizes:\n")
cat("Low Zn:  ", nrow(low_dat), "obs,", length(unique(low_dat$study_id)), "studies\n")
cat("High Zn: ", nrow(high_dat), "obs,", length(unique(high_dat$study_id)), "studies\n")
cat("All:     ", nrow(all_dat), "obs,", length(unique(all_dat$study_id)), "studies\n")

if(has_plant_part) {
  cat("\nPlant part distribution:\n")
  print(table(all_data$plant_part, all_data$Zn_Level, useNA = "ifany"))
}

# ============================================
# 6. Model fitting
# ============================================
cat("\n========== Fitting Models ==========\n")

fit_mod <- function(dat, mods_formula, name) {
  cat("\n>>> Fitting", name, "...\n")
  tryCatch({
    mod <- rma.mv(
      yi = lnRR_biomass,
      V = variance_biomass,
      mods = mods_formula,
      random = ~ 1 | study_id / obs_id,
      data = dat,
      method = "REML",
      control = list(optimizer = "bobyqa", maxiter = 10000, stepadj = 0.5)
    )
    cat("✅", name, "fitted\n")
    return(mod)
  }, error = function(e) {
    message("❌ ", name, " fitting failed: ", e$message)
    return(NULL)
  })
}

mod_low  <- fit_mod(low_dat, ~ lnRR_Zn + lnRR_P, "Low Zn")
mod_high <- fit_mod(high_dat, ~ lnRR_Zn + lnRR_P, "High Zn")
mod_all  <- fit_mod(all_dat, ~ lnRR_Zn + lnRR_P, "All")
mod_inter <- fit_mod(reg_data, 
                     ~ lnRR_Zn + lnRR_P + Zn_Level + lnRR_Zn:Zn_Level + lnRR_P:Zn_Level, 
                     "Interaction")

if(any(sapply(list(mod_low, mod_high, mod_all, mod_inter), is.null))) {
  warning("Some models failed to fit. Check output above.")
}

cat("\n========== Convergence Status ==========\n")
cat("Low Zn:    ", ifelse(!is.null(mod_low) && isTRUE(mod_low$converged), "✅", "❌"), "\n")
cat("High Zn:   ", ifelse(!is.null(mod_high) && isTRUE(mod_high$converged), "✅", "❌"), "\n")
cat("All:       ", ifelse(!is.null(mod_all) && isTRUE(mod_all$converged), "✅", "❌"), "\n")
cat("Interaction:", ifelse(!is.null(mod_inter) && isTRUE(mod_inter$converged), "✅", "❌"), "\n")

# ============================================
# 7. Result extraction and summary
# ============================================
cat("\nCalculating heterogeneity...\n")
hetero_low_fixed  <- report_heterogeneity_fixed(mod_low, "Low Zn")
hetero_high_fixed <- report_heterogeneity_fixed(mod_high, "High Zn")
hetero_all_fixed  <- report_heterogeneity_fixed(mod_all, "All")
hetero_fixed <- bind_rows(hetero_low_fixed, hetero_high_fixed, hetero_all_fixed)
cat("\n【Heterogeneity Results】\n")
print(hetero_fixed, digits = 4)

cat("\n========== Extracting Coefficients ==========\n")
coef_table <- bind_rows(
  extract_coef(mod_low, "Low Zn", apply_fdr = TRUE),
  extract_coef(mod_high, "High Zn", apply_fdr = TRUE),
  extract_coef(mod_all, "All", apply_fdr = TRUE)
)
cat("\n【Subgroup Model Coefficients - FDR adjusted】\n")
print(coef_table, digits = 4)

coef_table_inter <- extract_coef(mod_inter, "Interaction Model", apply_fdr = FALSE)
cat("\n【Interaction Model Coefficients - Unadjusted】\n")
print(coef_table_inter, digits = 4)

sample_table <- data.frame(
  Model = c("Low Zn", "High Zn", "All"),
  n_obs = c(nrow(low_dat), nrow(high_dat), nrow(all_dat)),
  n_study = c(length(unique(low_dat$study_id)),
              length(unique(high_dat$study_id)),
              length(unique(all_dat$study_id))),
  avg_obs_per_study = c(round(nrow(low_dat)/length(unique(low_dat$study_id)), 1),
                        round(nrow(high_dat)/length(unique(high_dat$study_id)), 1),
                        round(nrow(all_dat)/length(unique(all_dat$study_id)), 1))
)
cat("\n【Sample Info】\n")
print(sample_table)

summary_display <- coef_table %>%
  filter(Term != "intrcpt") %>%
  mutate(
    Term = ifelse(Term == "lnRR_Zn", "β_Zn", "β_P"),
    `Estimate±SE` = paste0(round(Estimate, 3), "±", round(SE, 3)),
    `p` = paste0(round(p_fdr, 4), Significance)
  ) %>%
  select(Model, Term, `Estimate±SE`, `p`) %>%
  pivot_wider(names_from = Term, values_from = c(`Estimate±SE`, `p`)) %>%
  rename(`β_Zn` = `Estimate±SE_β_Zn`, `p_Zn` = p_β_Zn,
         `β_P` = `Estimate±SE_β_P`, `p_P` = p_β_P)

summary_display_updated <- summary_display %>%
  left_join(hetero_fixed %>% select(Model, I2_study, I2_obs, I2_total), by = "Model") %>%
  mutate(
    `I²_study` = paste0(I2_study, "%"),
    `I²_obs` = paste0(I2_obs, "%"),
    `I²_total` = paste0(I2_total, "%")
  ) %>%
  left_join(sample_table, by = "Model") %>%
  mutate(`n (obs/study)` = paste0(n_obs, "/", n_study)) %>%
  select(Model, `n (obs/study)`, `β_Zn`, `p_Zn`, `β_P`, `p_P`, `I²_study`, `I²_obs`, `I²_total`)

cat("\n【Summary Table】\n")
print(summary_display_updated)

# Export tables
write.csv(coef_table, file.path(FIG_OUT_PATH, "B_coefficients_table.csv"), row.names = FALSE)
write.csv(coef_table_inter, file.path(FIG_OUT_PATH, "B_coefficients_interaction_table.csv"), row.names = FALSE)
write.csv(hetero_fixed, file.path(FIG_OUT_PATH, "B_heterogeneity_table.csv"), row.names = FALSE)
write.csv(sample_table, file.path(FIG_OUT_PATH, "B_sample_table.csv"), row.names = FALSE)
write.csv(summary_display_updated, file.path(FIG_OUT_PATH, "summary_display_fixed.csv"), row.names = FALSE)

# Save model objects (new)
saveRDS(list(low = mod_low, high = mod_high, all = mod_all, inter = mod_inter), 
        file.path(FIG_OUT_PATH, "meta_models.rds"))
cat("\n✅ Tables and models saved\n")

# ============================================
# 8. Plotting and export
# ============================================
cat("\n========== Plotting ==========\n")

# Figure 1: Scatter plots
cat("\n--- Figure 1: Scatter plots ---\n")
p1A <- create_scatter_with_stats(reg_data, "lnRR_Zn", "lnRR_biomass", 
                                 "lnRR_Zn", "lnRR_biomass", 
                                 "A. lnRR_Zn vs lnRR_biomass")
p1B <- create_scatter_with_stats(reg_data, "lnRR_P", "lnRR_biomass", 
                                 "lnRR_P", "lnRR_biomass", 
                                 "B. lnRR_P vs lnRR_biomass")
p1C <- create_scatter_with_stats(reg_data, "lnRR_Zn", "lnRR_P", 
                                 "lnRR_Zn", "lnRR_P", 
                                 "C. lnRR_Zn vs lnRR_P")

Fig1 <- (p1A | p1B | p1C) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

safe_save(Fig1, file.path(FIG_OUT_PATH, "B_Figure1_Scatter_with_Weight.pdf"), w = 15, h = 6)
safe_save(Fig1, file.path(FIG_OUT_PATH, "B_Figure1_Scatter_with_Weight.png"), w = 15, h = 6)

# Figure 2d: Subgroup coefficients
cat("\n--- Figure 2d: Subgroup coefficients ---\n")
coef_plot_data <- coef_table %>%
  filter(Term != "intrcpt") %>%
  mutate(
    Term_label = factor(ifelse(Term == "lnRR_Zn", "lnRR_Zn", "lnRR_P"),
                        levels = c("lnRR_Zn", "lnRR_P")),
    Model = factor(Model, levels = c("Low Zn", "High Zn", "All")),
    Sig_label = Significance
  )

ymax_B <- max(abs(c(coef_plot_data$Estimate + coef_plot_data$SE,
                    coef_plot_data$Estimate - coef_plot_data$SE)), na.rm = TRUE) * 1.4
if(ymax_B == 0 || is.na(ymax_B) || is.infinite(ymax_B)) ymax_B <- 1

p2d <- ggplot(coef_plot_data, aes(x = Term_label, y = Estimate, color = Model, shape = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_pointrange(aes(ymin = Estimate - SE, ymax = Estimate + SE),
                  position = position_dodge(width = 0.5),
                  linewidth = 1, size = 0.8, fatten = 3, stroke = 1.2) +
  geom_text(aes(label = Sig_label, 
                y = Estimate + sign(Estimate) * 0.1 * ymax_B),
            position = position_dodge(width = 0.5),
            size = 5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = GROUP_COLORS, name = "Zn Level") +
  scale_shape_manual(values = GROUP_SHAPES, name = "Zn Level") +
  labs(title = "D. Subgroup meta-regression coefficients (FDR-adjusted)",
       x = "Predictor", 
       y = expression("Regression coefficient (" * beta * ") ± SE"),
       caption = paste0("n = ", nrow(reg_data), ", k = ", length(unique(reg_data$study_id)), 
                        " | FDR-adjusted p-values: *p<0.05, **p<0.01, ***p<0.001")) +
  coord_cartesian(ylim = c(-ymax_B, ymax_B)) +
  theme_meta() +
  theme(panel.grid.major.x = element_blank(),
        legend.box = "horizontal")

safe_save(p2d, file.path(FIG_OUT_PATH, "B_Figure2d_Subgroup_Coefficients.pdf"), w = 8, h = 5.5)
safe_save(p2d, file.path(FIG_OUT_PATH, "B_Figure2d_Subgroup_Coefficients.png"), w = 8, h = 5.5)

# Figure 2e: Interaction coefficients
cat("\n--- Figure 2e: Interaction coefficients ---\n")
coef_plot_data_inter <- coef_table_inter %>%
  filter(Term != "intrcpt") %>%
  mutate(
    Term_display = case_when(
      Term == "lnRR_Zn" ~ "lnRR_Zn",
      Term == "lnRR_P" ~ "lnRR_P",
      Term == "Zn_LevelHigh Zn" ~ "Zn Level\n(High vs Low)",
      Term == "lnRR_Zn:Zn_LevelHigh Zn" ~ "lnRR_Zn ×\nZn Level",
      Term == "lnRR_P:Zn_LevelHigh Zn" ~ "lnRR_P ×\nZn Level",
      TRUE ~ Term
    ),
    Category = ifelse(str_detect(Term, ":"), "Interaction", "Main effect"),
    Sig_label = Significance,
    Term_display = fct_relevel(Term_display, "lnRR_Zn", "lnRR_P", "Zn Level\n(High vs Low)", 
                               "lnRR_Zn ×\nZn Level", "lnRR_P ×\nZn Level")
  )

ymax_inter <- max(abs(c(coef_plot_data_inter$Estimate + coef_plot_data_inter$SE,
                        coef_plot_data_inter$Estimate - coef_plot_data_inter$SE)), na.rm = TRUE) * 1.4
if(ymax_inter == 0 || is.na(ymax_inter) || is.infinite(ymax_inter)) ymax_inter <- 1

p2e <- ggplot(coef_plot_data_inter, aes(x = Term_display, y = Estimate, color = Category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_pointrange(aes(ymin = Estimate - SE, ymax = Estimate + SE),
                  linewidth = 0.8, size = 0.6, fatten = 2.5) +
  geom_text(aes(label = Sig_label, y = Estimate + sign(Estimate) * 0.1 * ymax_inter),
            size = 5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("Main effect" = "#2c3e50", "Interaction" = "#e74c3c")) +
  labs(title = "E. Interaction model coefficients (unadjusted p-values)",
       x = "Predictor", y = expression("Regression coefficient (" * beta * ") ± SE"),
       color = "Effect type",
       caption = paste0("n = ", nrow(reg_data), ", k = ", length(unique(reg_data$study_id)),
                        " | Unadjusted p-values: *p<0.05, **p<0.01, ***p<0.001")) +
  coord_cartesian(ylim = c(-ymax_inter, ymax_inter)) +
  theme_meta() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 9.5))

safe_save(p2e, file.path(FIG_OUT_PATH, "B_Figure2e_Interaction_Coefficients.pdf"), w = 10, h = 5.5)
safe_save(p2e, file.path(FIG_OUT_PATH, "B_Figure2e_Interaction_Coefficients.png"), w = 10, h = 5.5)

# Figure 2f: Mean effects
cat("\n--- Figure 2f: Mean effects ---\n")
p2f <- NULL

bio_all <- all_data %>% 
  filter(!is.na(lnRR_biomass), !is.na(variance_biomass)) %>% 
  rename(lnRR = lnRR_biomass, variance = variance_biomass)
zn_all <- all_data %>% 
  filter(!is.na(lnRR_Zn), !is.na(variance_Zn)) %>% 
  rename(lnRR = lnRR_Zn, variance = variance_Zn)
p_all <- all_data %>% 
  filter(!is.na(lnRR_P), !is.na(variance_P)) %>% 
  rename(lnRR = lnRR_P, variance = variance_P)

mean_res <- bind_rows(
  calc_mean_effect(bio_all %>% filter(Zn_Level == "Low Zn"), "Biomass", "Low Zn"),
  calc_mean_effect(bio_all %>% filter(Zn_Level == "High Zn"), "Biomass", "High Zn"),
  calc_mean_effect(zn_all %>% filter(Zn_Level == "Low Zn"), "Zn", "Low Zn"),
  calc_mean_effect(zn_all %>% filter(Zn_Level == "High Zn"), "Zn", "High Zn"),
  calc_mean_effect(p_all %>% filter(Zn_Level == "Low Zn"), "P", "Low Zn"),
  calc_mean_effect(p_all %>% filter(Zn_Level == "High Zn"), "P", "High Zn")
)

if(!is.null(mean_res) && nrow(mean_res) > 0){
  mean_res <- mean_res %>% filter(!is.na(p))
  if(nrow(mean_res) > 0){
    mean_res$p_fdr <- fdr_correct(mean_res$p)
    mean_res$Sig <- sig_star(mean_res$p_fdr)
    mean_res$Group <- factor(mean_res$Group, levels = c("Low Zn", "High Zn"))
    
    y_max <- max(mean_res$ci_ub, na.rm = TRUE)
    y_min <- min(mean_res$ci_lb, na.rm = TRUE)
    y_range <- diff(c(y_min, y_max))
    
    p2f <- ggplot(mean_res, aes(x = Response, y = mean_lnRR, fill = Group)) +
      geom_col(position = position_dodge(0.7), width = 0.6, color = "gray30", linewidth = 0.3) +
      geom_errorbar(aes(ymin = ci_lb, ymax = ci_ub), 
                    width = 0.2, position = position_dodge(0.7), linewidth = 0.7) +
      geom_text(aes(label = Sig, 
                    y = ifelse(mean_lnRR > 0, ci_ub + 0.05 * y_range, ci_lb - 0.05 * y_range)), 
                position = position_dodge(0.7), size = 5, fontface = "bold") +
      geom_text(aes(label = paste0("n=", n)), 
                position = position_dodge(0.7), 
                vjust = 1.8, size = 3, color = "gray40", fontface = "italic") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
      scale_fill_manual(values = GROUP_COLORS) +
      labs(title = "F. Overall mean effect sizes (FDR-adjusted)", 
           x = "Variable", y = "Mean lnRR ± 95% CI",
           caption = "FDR-adjusted p-values: *p<0.05, **p<0.01, ***p<0.001") +
      theme_meta() +
      theme(panel.grid.major.x = element_blank())
    
    safe_save(p2f, file.path(FIG_OUT_PATH, "B_Figure2f_MeanEffect.pdf"), w = 10, h = 5.5)
    safe_save(p2f, file.path(FIG_OUT_PATH, "B_Figure2f_MeanEffect.png"), w = 10, h = 5.5)
  }
}

# Figure 3: Summary table
cat("\n--- Figure 3: Summary table ---\n")
p3 <- create_table_plot(summary_display_updated, "Summary of model parameters", font_size = 4.2)
safe_save(p3, file.path(FIG_OUT_PATH, "B_Figure3_SummaryTable.pdf"), w = 10, h = 4)
safe_save(p3, file.path(FIG_OUT_PATH, "B_Figure3_SummaryTable.png"), w = 10, h = 4)

# Figure 4: Plant part comparison
Figure4 <- NULL
if(has_plant_part) {
  Figure4 <- create_figure4_optimized(all_data)
  if(!is.null(Figure4)) {
    safe_save(Figure4, file.path(FIG_OUT_PATH, "B_Figure4_Plant_Part_Comparison_Optimized.pdf"), w = 14, h = 6)
    safe_save(Figure4, file.path(FIG_OUT_PATH, "B_Figure4_Plant_Part_Comparison_Optimized.png"), w = 14, h = 6)
    cat("\n✅ Figure 4 saved\n")
  }
} else {
  cat("\n⚠️ Skipping Figure 4: missing plant_part column\n")
}

# ============================================
# 9. Combine overall figure and final output
# ============================================
cat("\n--- Combining final figure ---\n")

if(!is.null(p2f) && !is.null(Figure4)) {
  Fig_combined <- (Fig1) / (p2d | p2e) / (p2f) / (p3) / (Figure4) +
    plot_layout(heights = c(1, 1, 0.9, 0.5, 1.2)) +
    plot_annotation(
      title = "Zinc and Phosphorus Effects on Plant Biomass: A Meta-Analysis",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
  safe_save(Fig_combined, file.path(FIG_OUT_PATH, "B_Figure_Combined_with_Part.pdf"), w = 14, h = 20)
  safe_save(Fig_combined, file.path(FIG_OUT_PATH, "B_Figure_Combined_with_Part.png"), w = 14, h = 20)
} else if(!is.null(p2f)) {
  Fig_combined <- (Fig1) / (p2d | p2e) / (p2f) / (p3) +
    plot_layout(heights = c(1, 1, 0.9, 0.5)) +
    plot_annotation(
      title = "Zinc and Phosphorus Effects on Plant Biomass: A Meta-Analysis",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
  safe_save(Fig_combined, file.path(FIG_OUT_PATH, "B_Figure_Combined.pdf"), w = 14, h = 16)
  safe_save(Fig_combined, file.path(FIG_OUT_PATH, "B_Figure_Combined.png"), w = 14, h = 16)
} else {
  Fig_combined <- (Fig1) / (p2d | p2e) / (p3) +
    plot_layout(heights = c(1, 1, 0.5)) +
    plot_annotation(
      title = "Zinc and Phosphorus Effects on Plant Biomass: A Meta-Analysis",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
  safe_save(Fig_combined, file.path(FIG_OUT_PATH, "B_Figure_Combined.pdf"), w = 14, h = 12)
  safe_save(Fig_combined, file.path(FIG_OUT_PATH, "B_Figure_Combined.png"), w = 14, h = 12)
}

# ==================================================================
# Results summary output
# ==================================================================
cat("\n\n========== Final Results Summary ==========\n")
cat("\n【1】Coefficient estimates (FDR-adjusted)\n")
print(coef_table, digits = 4)
cat("\n【2】Interaction model coefficients (unadjusted)\n")
print(coef_table_inter, digits = 4)
cat("\n【3】Heterogeneity\n")
print(hetero_fixed, digits = 4)
cat("\n【4】Sample info\n")
print(sample_table)
cat("\n【5】Summary table\n")
print(summary_display_updated)

if(has_plant_part) {
  cat("\n【6】Plant part info\n")
  part_counts <- all_data %>%
    filter(!is.na(plant_part)) %>%
    group_by(Zn_Level, plant_part) %>%
    summarise(n = n(), .groups = "drop")
  print(part_counts)
}

# Extract key results
get_beta <- function(df, m, term){
  row <- df %>% filter(Model == m, Term == term)
  if(nrow(row) == 0) return(list(b = NA, p = "NA"))
  b_val <- round(row$Estimate, 3)
  pfdr <- row$p_fdr
  p_text <- if(pfdr < 0.001) "<0.001" else as.character(round(pfdr, 3))
  return(list(b = b_val, p = p_text))
}

p_low_zn  <- get_beta(coef_table, "Low Zn", "lnRR_Zn")
p_high_zn <- get_beta(coef_table, "High Zn", "lnRR_Zn")
p_all_zn  <- get_beta(coef_table, "All", "lnRR_Zn")
p_low_p   <- get_beta(coef_table, "Low Zn", "lnRR_P")
p_high_p  <- get_beta(coef_table, "High Zn", "lnRR_P")
p_all_p   <- get_beta(coef_table, "All", "lnRR_P")

cat("\n\n========== Key Findings ==========\n")
cat("\n📊 Phosphorus (P) effect:\n")
cat("   - Low Zn:  β =", p_low_p$b,  ", FDR-p =", p_low_p$p, "\n")
cat("   - High Zn: β =", p_high_p$b, ", FDR-p =", p_high_p$p, "\n")
cat("   - All:     β =", p_all_p$b,  ", FDR-p =", p_all_p$p, "\n")

cat("\n📊 Zinc (Zn) effect:\n")
cat("   - Low Zn:  β =", p_low_zn$b,  ", FDR-p =", p_low_zn$p, "\n")
cat("   - High Zn: β =", p_high_zn$b, ", FDR-p =", p_high_zn$p, "\n")
cat("   - All:     β =", p_all_zn$b,  ", FDR-p =", p_all_zn$p, "\n")

cat("\n📊 Hierarchical heterogeneity I²:\n")
print(hetero_fixed[, c("Model", "I2_study", "I2_obs", "I2_total")])

diagnose_subgroup(mod_low, low_dat, "Low Zn")
diagnose_subgroup(mod_high, high_dat, "High Zn")
diagnose_subgroup(mod_all, all_dat, "All")

cat("\n\n✅ Analysis complete. Output directory:", normalizePath(FIG_OUT_PATH), "\n")

