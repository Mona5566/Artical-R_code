# ============================================
# Meta regression combined panel plot (fully automatic combination + results table)
# Supports Total_Zn and soil_pH
# Output: combined panel PDF + results summary table CSV/PDF
# ============================================

library(tidyverse)
library(ggplot2)
library(metafor)
library(patchwork)
library(cowplot)
library(gridExtra)   # For table PDF

# ============================================
# Configuration
# ============================================
PROCESSED_DIR <- "data/processed"
FIGURES_DIR   <- "figures"

ZN_LEVEL_COLORS <- c("Low Zn" = "#2E86AB", "High Zn" = "#A23B72")

# ============================================
# Theme
# ============================================
theme_metareg <- function() {
  theme_bw() %+replace%
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.title       = element_text(face = "bold", size = 8),
      axis.text        = element_text(color = "black", size = 7),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "gray30", fill = NA, linewidth = 0.6),
      legend.position  = "none",
      legend.key       = element_blank(),
      legend.title     = element_text(size = 7, face = "bold"),
      legend.text      = element_text(size = 6),
      complete         = TRUE
    )
}

# ============================================
# Get category
# ============================================
get_category <- function(x) {
  if(grepl("Biomass", x, ignore.case = TRUE)) return("Biomass")
  if(grepl("Zn_|Zn$|Zinc|Zn concentration", x, ignore.case = TRUE)) return("Zn_concentration")
  if(grepl("BCF", x, ignore.case = TRUE)) return("BCF")
  if(grepl("P_|Phosphorus|P$|Phos|Total_P|P content", x, ignore.case = TRUE)) return("Phosphorus")
  return("Other")
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
# Single panel (returns plot + result table row)
# ============================================
plot_metareg_single <- function(dat_all, dat_low, dat_high,
                                resp, moderator = "Total_Zn",
                                tag = NULL) {
  
  empty_result <- data.frame(
    Moderator = moderator, Response = resp, Category = get_category(resp),
    n_All = NA, n_Low = NA, n_High = NA,
    slope_All = NA, p_All = NA, sig_All = NA,
    slope_Low = NA, p_Low = NA, sig_Low = NA,
    slope_High = NA, p_High = NA, sig_High = NA,
    Plotted = FALSE,
    stringsAsFactors = FALSE
  )
  
  if(is.null(dat_all) || !resp %in% names(dat_all)) return(list(plot = NULL, result = empty_result))
  
  df_all  <- prepare_moderator_data(dat_all[[resp]],  moderator)
  df_low  <- if(!is.null(dat_low)  && resp %in% names(dat_low))  prepare_moderator_data(dat_low[[resp]],  moderator) else NULL
  df_high <- if(!is.null(dat_high) && resp %in% names(dat_high)) prepare_moderator_data(dat_high[[resp]], moderator) else NULL
  
  if(is.null(df_all)) return(list(plot = NULL, result = empty_result))
  
  model_all  <- fit_metareg_model(df_all)
  model_low  <- if(!is.null(df_low))  fit_metareg_model(df_low)  else NULL
  model_high <- if(!is.null(df_high)) fit_metareg_model(df_high) else NULL
  
  if(is.null(model_all)) return(list(plot = NULL, result = empty_result))
  
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
  
  # Scatter
  low_points <- if(!is.null(df_low) && nrow(df_low) > 0) {
    data.frame(moderator_trans = df_low$moderator_trans,
               lnRR = df_low$lnRR,
               precision = 1 / sqrt(df_low$variance),
               group = "Low Zn")
  } else NULL
  
  high_points <- if(!is.null(df_high) && nrow(df_high) > 0) {
    data.frame(moderator_trans = df_high$moderator_trans,
               lnRR = df_high$lnRR,
               precision = 1 / sqrt(df_high$variance),
               group = "High Zn")
  } else NULL
  
  all_points <- data.frame(moderator_trans = df_all$moderator_trans,
                           lnRR = df_all$lnRR,
                           precision = 1 / sqrt(df_all$variance))
  
  color_points <- rbind(low_points, high_points)
  if(is.null(color_points) || nrow(color_points) == 0) {
    return(list(plot = NULL, result = empty_result))
  }
  color_points$group <- factor(color_points$group, levels = c("Low Zn", "High Zn"))
  
  # Precision breaks
  all_precision <- c(all_points$precision, color_points$precision)
  precision_breaks <- quantile(all_precision, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  if(length(unique(round(precision_breaks, 1))) < 2) {
    precision_breaks <- c(min(all_precision, na.rm = TRUE), max(all_precision, na.rm = TRUE))
  }
  
  # Statistical labels
  stats_text <- ""
  if(!is.null(model_low)) {
    stats_text <- paste0(stats_text, "Low: ", sprintf("%.3f", model_low$b[2]),
                         ifelse(model_low$pval[2] < 0.05, "*", ""), "  ")
  }
  if(!is.null(model_high)) {
    stats_text <- paste0(stats_text, "High: ", sprintf("%.3f", model_high$b[2]),
                         ifelse(model_high$pval[2] < 0.05, "*", ""), "  ")
  }
  if(!is.null(model_all)) {
    stats_text <- paste0(stats_text, "All: ", sprintf("%.3f", model_all$b[2]),
                         ifelse(model_all$pval[2] < 0.05, "*", ""))
  }
  
  x_label <- get_moderator_label(moderator)
  
  p <- ggplot() +
    geom_ribbon(data = pred_all_df, aes(x = x, ymin = ci_lb, ymax = ci_ub),
                fill = "gray50", alpha = 0.15) +
    geom_line(data = pred_all_df, aes(x = x, y = y),
              color = "gray40", linewidth = 0.7, linetype = "dashed") +
    {if(!is.null(pred_low_df)) geom_ribbon(data = pred_low_df, aes(x = x, ymin = ci_lb, ymax = ci_ub),
                                           fill = ZN_LEVEL_COLORS["Low Zn"], alpha = 0.15)} +
    {if(!is.null(pred_low_df)) geom_line(data = pred_low_df, aes(x = x, y = y),
                                         color = ZN_LEVEL_COLORS["Low Zn"], linewidth = 1.0)} +
    {if(!is.null(pred_high_df)) geom_ribbon(data = pred_high_df, aes(x = x, ymin = ci_lb, ymax = ci_ub),
                                            fill = ZN_LEVEL_COLORS["High Zn"], alpha = 0.15)} +
    {if(!is.null(pred_high_df)) geom_line(data = pred_high_df, aes(x = x, y = y),
                                          color = ZN_LEVEL_COLORS["High Zn"], linewidth = 1.0)} +
    geom_point(data = all_points,
               aes(x = moderator_trans, y = lnRR, size = precision),
               color = "gray60", alpha = 0.2, shape = 1) +
    geom_point(data = color_points,
               aes(x = moderator_trans, y = lnRR, size = precision, color = group),
               alpha = 0.85) +
    scale_color_manual(values = ZN_LEVEL_COLORS) +
    scale_size_continuous(range = c(1.2, 5), breaks = round(precision_breaks, 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
    annotate("text", x = -Inf, y = Inf,
             label = paste("Slopes:\n", stats_text),
             hjust = -0.05, vjust = 1.3, size = 2.2, fontface = "italic", color = "gray30") +
    annotate("text", x = -Inf, y = Inf,
             label = paste0("n_All=", nrow(df_all),
                            ifelse(!is.null(df_low), paste0(", n_Low=", nrow(df_low)), ""),
                            ifelse(!is.null(df_high), paste0(", n_High=", nrow(df_high)), "")),
             hjust = -0.05, vjust = 0.1, size = 2.0, color = "gray40") +
    labs(title = resp, x = x_label, y = "lnRR") +
    theme_metareg()
  
  if(!is.null(tag)) {
    p <- p + labs(tag = tag) +
      theme(plot.tag = element_text(face = "bold", size = 11),
            plot.tag.position = c(0, 1))
  }
  
  # Result table row
  res <- data.frame(
    Moderator = moderator,
    Response  = resp,
    Category  = get_category(resp),
    n_All     = nrow(df_all),
    n_Low     = ifelse(is.null(df_low), NA, nrow(df_low)),
    n_High    = ifelse(is.null(df_high), NA, nrow(df_high)),
    slope_All  = ifelse(is.null(model_all), NA, round(model_all$b[2], 4)),
    p_All      = ifelse(is.null(model_all), NA, round(model_all$pval[2], 4)),
    sig_All    = ifelse(is.null(model_all), NA, ifelse(model_all$pval[2] < 0.05, "*", "")),
    slope_Low  = ifelse(is.null(model_low), NA, round(model_low$b[2], 4)),
    p_Low      = ifelse(is.null(model_low), NA, round(model_low$pval[2], 4)),
    sig_Low    = ifelse(is.null(model_low), NA, ifelse(model_low$pval[2] < 0.05, "*", "")),
    slope_High = ifelse(is.null(model_high), NA, round(model_high$b[2], 4)),
    p_High     = ifelse(is.null(model_high), NA, round(model_high$pval[2], 4)),
    sig_High   = ifelse(is.null(model_high), NA, ifelse(model_high$pval[2] < 0.05, "*", "")),
    Plotted    = TRUE,
    stringsAsFactors = FALSE
  )
  
  return(list(plot = p, result = res))
}

# ============================================
# Fully automatic combination + results table
# ============================================
run_all_combined <- function(moderators = c("Total_Zn", "soil_pH"),
                             ncol = 4,
                             out_width = 16,
                             sub_height = 4) {
  
  cat("\n========================================\n")
  cat("  Fully automatic Meta regression combined panel + results table\n")
  cat("========================================\n\n")
  
  file_all  <- file.path(PROCESSED_DIR, "all_processed_data.rds")
  file_low  <- file.path(PROCESSED_DIR, "all_processed_data_low.rds")
  file_high <- file.path(PROCESSED_DIR, "all_processed_data_high.rds")
  
  if(!all(file.exists(c(file_all, file_low, file_high)))) {
    stop("Missing All / Low / High data files, please check data/processed directory.")
  }
  
  data_all  <- readRDS(file_all)
  data_low  <- readRDS(file_low)
  data_high <- readRDS(file_high)
  
  # Automatically get all response variables
  exclude_vars <- c("study_id", "variance", "lnRR", "Total_Zn", "soil_pH",
                    "moderator_transformed", "log_zn", "moderator_trans",
                    "moderator_orig", "obs_id")
  responses <- setdiff(names(data_all), exclude_vars)
  
  cat(sprintf("Found %d candidate response variables\n", length(responses)))
  
  all_results <- list()
  all_tables  <- list()
  
  for(mod in moderators) {
    cat(sprintf("\n========== Moderator: %s ==========\n", mod))
    
    plots_list <- list()
    table_rows <- list()
    
    for(i in seq_along(responses)) {
      resp <- responses[i]
      cat(sprintf("Processing: %s | %s\n", resp, mod))
      
      out <- plot_metareg_single(
        dat_all  = data_all,
        dat_low  = data_low,
        dat_high = data_high,
        resp     = resp,
        moderator = mod,
        tag      = letters[i]
      )
      
      table_rows[[resp]] <- out$result
      
      if(!is.null(out$plot)) {
        plots_list[[resp]] <- out$plot
        cat("  ✅ Done\n")
      } else {
        cat("  ⚠ Skipped\n")
      }
    }
    
    # Merge result table
    mod_table <- bind_rows(table_rows) %>%
      arrange(Category, Response)
    all_tables[[mod]] <- mod_table
    
    # Save result table CSV
    out_dir <- file.path(FIGURES_DIR, "metareg_plots", paste0("panels_combined_", mod))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    write.csv(mod_table,
              file.path(out_dir, paste0("metareg_results_", mod, ".csv")),
              row.names = FALSE)
    
    # Save result table PDF (table format)
    pdf(file.path(out_dir, paste0("metareg_results_", mod, ".pdf")),
        width = 16, height = max(4, 0.35 * nrow(mod_table) + 2))
    grid.extra <- tableGrob(
      mod_table %>%
        select(Moderator, Response, Category, n_All, n_Low, n_High,
               slope_All, p_All, sig_All,
               slope_Low, p_Low, sig_Low,
               slope_High, p_High, sig_High),
      rows = NULL,
      theme = ttheme_minimal(base_size = 7)
    )
    grid.arrange(grid.extra, top = paste("Meta-regression results:", mod))
    dev.off()
    
    cat(sprintf("  💾 Result table: metareg_results_%s.csv / .pdf\n", mod))
    
    # Combined panel
    if(length(plots_list) > 0) {
      nrow_plots <- ceiling(length(plots_list) / ncol)
      panel <- wrap_plots(plots_list, ncol = ncol, nrow = nrow_plots) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom",
              legend.title = element_text(size = 8, face = "bold"),
              legend.text  = element_text(size = 7))
      
      panel_file <- file.path(out_dir, paste0("panel_combined_", mod, ".pdf"))
      ggsave(panel_file, panel,
             width = out_width,
             height = sub_height * nrow_plots,
             device = "pdf")
      cat(sprintf("  💾 Combined plot: %s\n", panel_file))
      
      all_results[[mod]] <- panel
    } else {
      cat("  ⚠ No plots to combine\n")
    }
  }
  
  # Merge result tables for all moderators
  final_table <- bind_rows(all_tables)
  write.csv(final_table,
            file.path(FIGURES_DIR, "metareg_plots", "metareg_results_ALL.csv"),
            row.names = FALSE)
  
  cat("\n========================================\n")
  cat("✅ All complete!\n")
  cat("Result table: figures/metareg_plots/metareg_results_ALL.csv\n")
  cat("========================================\n")
  
  return(list(plots = all_results, table = final_table))
}

# ============================================
# Execute
# ============================================
results <- run_all_combined(
  moderators  = c("Total_Zn", "soil_pH"),
  ncol        = 4,
  out_width   = 16,
  sub_height  = 4
)