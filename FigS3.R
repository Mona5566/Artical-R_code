# ============================================
# Meta regression combined panel plot (all variables, 13 panels complete version)
# Supports Total_Zn and soil_pH moderators
# Low/High/All three regression lines + scatter + panel labels
# ============================================

library(tidyverse)
library(ggplot2)
library(metafor)
library(patchwork)
library(cowplot)

# ============================================
# Configuration
# ============================================
PROCESSED_DIR <- "data/processed"
FIGURES_DIR   <- "figures"

ZN_LEVEL_COLORS <- c("Low Zn" = "#2E86AB", "High Zn" = "#A23B72")

COLOR_MAP <- c(
  "Biomass"          = "#2E86AB",
  "Zn_concentration" = "#A23B72",
  "BCF"              = "#F18F01",
  "Phosphorus"       = "#4CAF50"
)

# ============================================
# Theme (adapted for compact multi-panel style)
# ============================================
theme_metareg <- function() {
  theme_bw() %+replace%
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.title      = element_text(face = "bold", size = 8),
      axis.text       = element_text(color = "black", size = 7),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border    = element_rect(color = "gray30", fill = NA, linewidth = 0.6),
      legend.position = "none",   # Combined plot has unified legend
      legend.key      = element_blank(),
      legend.title    = element_text(size = 8, face = "bold"),
      legend.text     = element_text(size = 7),
      complete        = TRUE
    )
}

# ============================================
# Moderator axis labels
# ============================================
get_moderator_label <- function(moderator) {
  if(moderator == "Total_Zn") return("Total Soil Zn (log10 mg kg⁻¹)")
  if(moderator == "soil_pH")  return("Soil pH")
  return(moderator)
}

# ============================================
# Preprocess data
# ============================================
prepare_moderator_data <- function(df, moderator) {
  if(!moderator %in% colnames(df)) return(NULL)
  
  df <- df %>%
    filter(!is.na(!!sym(moderator)), !is.na(lnRR),
           is.finite(lnRR), variance > 0)
  
  if(nrow(df) < 5) return(NULL)
  
  if(moderator == "Total_Zn") {
    df <- df %>% filter(!!sym(moderator) > 0)
    if(nrow(df) < 5) return(NULL)
    df$moderator_trans <- log10(df[[moderator]])
    df$moderator_orig  <- df[[moderator]]
  } else {
    df$moderator_trans <- df[[moderator]]
    df$moderator_orig  <- df[[moderator]]
  }
  
  if(length(unique(df$moderator_trans)) < 3) return(NULL)
  return(df)
}

# ============================================
# Fit model
# ============================================
fit_metareg_model <- function(df) {
  df$obs_id <- 1:nrow(df)
  
  model <- tryCatch({
    rma.mv(
      yi = lnRR, V = variance,
      mods = ~ moderator_trans,
      random = list(~ 1 | study_id, ~ 1 | obs_id),
      data = df, method = "REML"
    )
  }, error = function(e) {
    tryCatch({
      rma.mv(
        yi = lnRR, V = variance,
        mods = ~ moderator_trans,
        random = ~ 1 | study_id,
        data = df, method = "REML"
      )
    }, error = function(e2) NULL)
  })
  
  return(model)
}

# ============================================
# Single panel (no legend, for combination)
# ============================================
plot_metareg_single <- function(dat_all, dat_low, dat_high,
                                resp, moderator = "Total_Zn",
                                tag = NULL, global_breaks = NULL) {
  
  if(is.null(dat_all) || !resp %in% names(dat_all)) return(NULL)
  
  df_all  <- prepare_moderator_data(dat_all[[resp]],  moderator)
  df_low  <- if(!is.null(dat_low)  && resp %in% names(dat_low))  prepare_moderator_data(dat_low[[resp]],  moderator) else NULL
  df_high <- if(!is.null(dat_high) && resp %in% names(dat_high)) prepare_moderator_data(dat_high[[resp]], moderator) else NULL
  
  if(is.null(df_all)) return(NULL)
  
  model_all  <- fit_metareg_model(df_all)
  model_low  <- if(!is.null(df_low))  fit_metareg_model(df_low)  else NULL
  model_high <- if(!is.null(df_high)) fit_metareg_model(df_high) else NULL
  
  if(is.null(model_all)) return(NULL)
  
  # Prediction
  x_range <- seq(min(df_all$moderator_trans), max(df_all$moderator_trans), length.out = 100)
  pred_all <- predict(model_all, newmods = x_range, addx = TRUE)
  pred_all_df <- data.frame(x = x_range, y = pred_all$pred,
                            ci_lb = pred_all$ci.lb, ci_ub = pred_all$ci.ub)
  
  pred_low_df <- NULL
  if(!is.null(model_low) && !is.null(df_low)) {
    x_low <- seq(min(df_low$moderator_trans), max(df_low$moderator_trans), length.out = 100)
    pred_low <- predict(model_low, newmods = x_low, addx = TRUE)
    pred_low_df <- data.frame(x = x_low, y = pred_low$pred,
                              ci_lb = pred_low$ci.lb, ci_ub = pred_low$ci.ub)
  }
  
  pred_high_df <- NULL
  if(!is.null(model_high) && !is.null(df_high)) {
    x_high <- seq(min(df_high$moderator_trans), max(df_high$moderator_trans), length.out = 100)
    pred_high <- predict(model_high, newmods = x_high, addx = TRUE)
    pred_high_df <- data.frame(x = x_high, y = pred_high$pred,
                               ci_lb = pred_high$ci.lb, ci_ub = pred_high$ci.ub)
  }
  
  # Scatter data
  low_points <- if(!is.null(df_low) && nrow(df_low) > 0) {
    data.frame(moderator_trans = df_low$moderator_trans, lnRR = df_low$lnRR,
               precision = 1 / sqrt(df_low$variance), group = "Low Zn")
  } else NULL
  
  high_points <- if(!is.null(df_high) && nrow(df_high) > 0) {
    data.frame(moderator_trans = df_high$moderator_trans, lnRR = df_high$lnRR,
               precision = 1 / sqrt(df_high$variance), group = "High Zn")
  } else NULL
  
  all_points <- data.frame(moderator_trans = df_all$moderator_trans,
                           lnRR = df_all$lnRR,
                           precision = 1 / sqrt(df_all$variance))
  
  color_points <- rbind(low_points, high_points)
  if(is.null(color_points) || nrow(color_points) == 0) {
    color_points <- data.frame(moderator_trans = numeric(0), lnRR = numeric(0), 
                               precision = numeric(0), group = character(0))
  } else {
    color_points$group <- factor(color_points$group, levels = c("Low Zn", "High Zn"))
  }
  
  # Statistical labels
  stats_text <- ""
  if(!is.null(model_low)) {
    stats_text <- paste0(stats_text, "Low: ", sprintf("%.3f", model_low$b[2]), ifelse(model_low$pval[2] < 0.05, "*", ""), "  ")
  }
  if(!is.null(model_high)) {
    stats_text <- paste0(stats_text, "High: ", sprintf("%.3f", model_high$b[2]), ifelse(model_high$pval[2] < 0.05, "*", ""), "  ")
  }
  if(!is.null(model_all)) {
    stats_text <- paste0(stats_text, "All: ", sprintf("%.3f", model_all$b[2]), ifelse(model_all$pval[2] < 0.05, "*", ""))
  }
  
  x_label <- get_moderator_label(moderator)
  
  # Plot (key modification: color and size share same legend combination)
  p <- ggplot() +
    geom_ribbon(data = pred_all_df, aes(x = x, ymin = ci_lb, ymax = ci_ub), fill = "gray50", alpha = 0.15) +
    geom_line(data = pred_all_df, aes(x = x, y = y), color = "gray40", linewidth = 0.7, linetype = "dashed") +
    {if(!is.null(pred_low_df)) geom_ribbon(data = pred_low_df, aes(x = x, ymin = ci_lb, ymax = ci_ub), fill = ZN_LEVEL_COLORS["Low Zn"], alpha = 0.15)} +
    {if(!is.null(pred_low_df)) geom_line(data = pred_low_df, aes(x = x, y = y), color = ZN_LEVEL_COLORS["Low Zn"], linewidth = 1.0)} +
    {if(!is.null(pred_high_df)) geom_ribbon(data = pred_high_df, aes(x = x, ymin = ci_lb, ymax = ci_ub), fill = ZN_LEVEL_COLORS["High Zn"], alpha = 0.15)} +
    {if(!is.null(pred_high_df)) geom_line(data = pred_high_df, aes(x = x, y = y), color = ZN_LEVEL_COLORS["High Zn"], linewidth = 1.0)} +
    # All gray points: no legend
    geom_point(data = all_points, aes(x = moderator_trans, y = lnRR, size = precision), 
               color = "gray60", alpha = 0.2, shape = 1, show.legend = FALSE) +
    # Low/High colored points: put color and size in same aes mapping
    geom_point(data = color_points, aes(x = moderator_trans, y = lnRR, size = precision, color = group), 
               alpha = 0.85) +
    scale_color_manual(values = ZN_LEVEL_COLORS, name = "Zn Level") +
    scale_size_continuous(range = c(1.2, 5), breaks = global_breaks, name = "Precision") +
    # [Core] Force merge color and size legends, and make colored points slightly larger in legend
    guides(
      color = guide_legend(title = "Zn Level", override.aes = list(size = 3)),
      size  = guide_legend(title = "Precision")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
    annotate("text", x = -Inf, y = Inf, label = paste("Slopes:\n", stats_text), hjust = -0.05, vjust = 1.3, size = 2.2, fontface = "italic", color = "gray30") +
    annotate("text", x = -Inf, y = Inf, label = paste0("n_All=", nrow(df_all), ifelse(!is.null(df_low), paste0(", n_Low=", nrow(df_low)), ""), ifelse(!is.null(df_high), paste0(", n_High=", nrow(df_high)), "")), hjust = -0.05, vjust = 0.1, size = 2.0, color = "gray40") +
    labs(title = resp, x = x_label, y = "lnRR") +
    theme_metareg() +
    theme(legend.position = "none") # Hide within panel, rely on patchwork collection
  
  if(!is.null(tag)) {
    p <- p + labs(tag = tag) +
      theme(plot.tag = element_text(face = "bold", size = 11), plot.tag.position = c(0, 1))
  }
  
  return(p)
}

# ============================================
# ============================================
# Generate combined panel + perfect unified dual legend
# ============================================
run_panel_combined <- function(moderator = "Total_Zn",
                               response_vars = NULL,
                               ncol = 4,
                               out_width = 16,
                               out_height = 14) {
  
  cat("\n========================================\n")
  cat("  Generating Meta regression combined panel (perfect legend version)\n")
  cat("  Moderator:", moderator, "\n")
  cat("========================================\n\n")
  
  file_all  <- file.path(PROCESSED_DIR, "all_processed_data.rds")
  file_low  <- file.path(PROCESSED_DIR, "all_processed_data_low.rds")
  file_high <- file.path(PROCESSED_DIR, "all_processed_data_high.rds")
  
  data_all  <- readRDS(file_all)
  data_low  <- readRDS(file_low)
  data_high <- readRDS(file_high)
  
  # 1. Compute global Precision breaks
  cat("Computing global Precision breaks...\n")
  all_precisions <- c()
  for(resp in response_vars) {
    df_a <- if(resp %in% names(data_all)) prepare_moderator_data(data_all[[resp]], moderator) else NULL
    if(!is.null(df_a)) all_precisions <- c(all_precisions, 1 / sqrt(df_a$variance))
    df_l <- if(resp %in% names(data_low)) prepare_moderator_data(data_low[[resp]], moderator) else NULL
    if(!is.null(df_l)) all_precisions <- c(all_precisions, 1 / sqrt(df_l$variance))
    df_h <- if(resp %in% names(data_high)) prepare_moderator_data(data_high[[resp]], moderator) else NULL
    if(!is.null(df_h)) all_precisions <- c(all_precisions, 1 / sqrt(df_h$variance))
  }
  
  # Ensure extracted values are pure numbers without percentage names
  global_precision_breaks <- unname(quantile(all_precisions, probs = c(0.1, 0.5, 0.9), na.rm = TRUE))
  if(length(unique(round(global_precision_breaks, 1))) < 2) {
    global_precision_breaks <- c(min(all_precisions, na.rm = TRUE), max(all_precisions, na.rm = TRUE))
  }
  global_precision_breaks <- round(global_precision_breaks, 1)
  cat("Global Size breaks:", paste(global_precision_breaks, collapse = ", "), "\n\n")
  
  # 2. Draw panels
  plots_list <- list()
  for(i in seq_along(response_vars)) {
    resp <- response_vars[i]
    cat(sprintf("Processing: %s | %s\n", resp, moderator))
    p <- plot_metareg_single(
      dat_all = data_all, dat_low = data_low, dat_high = data_high,
      resp = resp, moderator = moderator, tag = letters[i], global_breaks = global_precision_breaks
    )
    if(!is.null(p)) {
      plots_list[[resp]] <- p
      cat("  ✅ Done\n")
    } else {
      cat("  ⚠ Skipped\n")
    }
  }
  
  # 3. Assemble panel (no collect, avoid duplicate legends)
  nrow_plots <- ceiling(length(plots_list) / ncol)
  panel_plots <- wrap_plots(plots_list, ncol = ncol, nrow = nrow_plots)
  
  # 4. Build a "perfect" dummy plot to extract legend
  # Dummy plot must have same data mapping to successfully extract color and size
  dummy_df <- data.frame(
    x = 1:6,
    y = 1:6,
    group = factor(rep(c("Low Zn", "High Zn"), each = 3), levels = c("Low Zn", "High Zn")),
    # Simulate same size structure as real plot so breaks display correctly
    precision = rep(global_precision_breaks, 2) 
  )
  
  p_dummy <- ggplot(dummy_df, aes(x = x, y = y, size = precision, color = group)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = ZN_LEVEL_COLORS, name = "Zn Level") +
    scale_size_continuous(range = c(1.2, 5), breaks = global_precision_breaks, name = "Precision") +
    guides(
      color = guide_legend(title = "Zn Level", override.aes = list(size = 4)),
      size  = guide_legend(title = "Precision")
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal", # Horizontal side-by-side
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8)
    )
  
  # Extract perfect legend from dummy plot
  unified_legend <- cowplot::get_legend(p_dummy)
  
  # 5. Final combination: panels on top + legend at bottom
  final_panel <- cowplot::plot_grid(
    panel_plots, 
    unified_legend, 
    ncol = 1, 
    rel_heights = c(1, 0.05) # Legend takes only 5% height
  )
  
  # Save
  out_dir <- file.path(FIGURES_DIR, "metareg_plots", paste0("panels_combined_", moderator))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0("panel_combined_", moderator, ".pdf"))
  
  ggsave(out_file, final_panel, width = out_width, height = out_height, device = "pdf")
  cat(sprintf("\n💾 Saved: %s\n", out_file))
  
  return(final_panel)
}

# ============================================
# Execute: integrate all 13 variables
# ============================================
# Extracted 13 response variables actually present in data frame
responses_all <- c(
  "Zn_shoot", "Zn_root", "Zn_plant", 
  "P_shoot", "P_root",
  "Biomass_shoot", "Biomass_root", "Biomass_fruit", "Biomass_leave", 
  "BCF_shoot", "BCF_root", "BCF_plant", 
  "Plant_height"
)

# 1. Total_Zn panel
p_total_zn <- run_panel_combined(
  moderator     = "Total_Zn",
  response_vars = responses_all,
  ncol          = 4,       # 4 columns: exactly 4 rows (4x4), height 14 to prevent squeezing
  out_width     = 16,
  out_height    = 14
)

# 2. soil_pH panel
p_ph <- run_panel_combined(
  moderator     = "soil_pH",
  response_vars = responses_all,
  ncol          = 4,       # Same 4-column layout
  out_width     = 16,
  out_height    = 14
)