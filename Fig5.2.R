# ============================================
# Fig5 redraw script: plot Meta regression fit by best model
# Standalone version (includes modeling functions, no need to source original script)
# ============================================

library(tidyverse)
library(metafor)
library(patchwork)

# ---------- 0. Path settings ----------
ROOT_DIR      <- "D:/2026-03_oncemore/TryRun22-40"
PROCESSED_DIR <- file.path(ROOT_DIR, "data", "processed")
FIGURES_DIR   <- file.path(ROOT_DIR, "results", "figures")
TABLES_DIR    <- file.path(ROOT_DIR, "results", "tables")

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------- 1. Read raw data ----------
data_file <- file.path(PROCESSED_DIR, "all_processed_data.rds")
if (!file.exists(data_file)) stop("❌ Data file not found: ", data_file)
processed_data <- readRDS(data_file)
cat("✅ Data loaded successfully:", data_file, "\n")

# ============================================
# 2. Core modeling functions (copied from original script)
# ============================================
fit_model <- function(formula_mods, data, use_random = TRUE) {
  if (use_random && length(unique(data$study_id)) >= 3) {
    tryCatch({
      return(rma.mv(
        yi = lnRR, V = variance,
        random = ~ 1 | study_id,
        mods = formula_mods,
        data = data,
        method = "REML",
        control = list(optimizer = "optim", optmethod = "Nelder-Mead", maxit = 2000)
      ))
    }, error = function(e) {
      tryCatch({
        return(rma(yi = lnRR, vi = variance, mods = formula_mods,
                   data = data, method = "REML"))
      }, error = function(e2) NULL)
    })
  } else {
    tryCatch({
      return(rma(yi = lnRR, vi = variance, mods = formula_mods,
                 data = data, method = "REML"))
    }, error = function(e) NULL)
  }
}

compare_metareg_models <- function(dat, moderator, response_name) {
  if (moderator == "Total_Zn") {
    if (any(dat[[moderator]] <= 0, na.rm = TRUE)) {
      dat$mod_transformed <- log10(dat[[moderator]] + 0.01)
    } else {
      dat$mod_transformed <- log10(dat[[moderator]])
    }
    x_label <- expression("Total soil Zn (mg/kg, log"[10] * ")")
  } else {
    dat$mod_transformed <- dat[[moderator]]
    x_label <- moderator
  }
  
  dat_clean <- dat %>%
    filter(!is.na(mod_transformed), !is.na(lnRR), is.finite(lnRR),
           variance > 0, is.finite(variance), !is.na(study_id)) %>%
    group_by(study_id) %>% filter(n() >= 2) %>% ungroup()
  
  if (nrow(dat_clean) < 10) return(NULL)
  if (length(unique(dat_clean$mod_transformed)) < 3) return(NULL)
  
  n_studies <- length(unique(dat_clean$study_id))
  
  # Linear
  model_linear <- fit_model(~ mod_transformed, dat_clean, use_random = n_studies >= 3)
  if (is.null(model_linear)) return(NULL)
  aic_linear <- tryCatch(AIC(model_linear), error = function(e) NA)
  p_linear   <- tryCatch(model_linear$pval[2], error = function(e) NA)
  linear_valid <- !is.na(p_linear) && p_linear < 0.05
  
  # Quadratic
  model_quad <- NULL; aic_quad <- NA; p_quad <- NA; quad_valid <- FALSE
  vertex <- NA; vertex_type <- NA
  if (nrow(dat_clean) >= 20) {
    x_mean <- mean(dat_clean$mod_transformed, na.rm = TRUE)
    dat_clean$mod_centered <- dat_clean$mod_transformed - x_mean
    dat_clean$mod_sq <- dat_clean$mod_centered^2
    cor_check <- tryCatch(cor(dat_clean$mod_centered, dat_clean$mod_sq, use = "complete.obs"),
                          error = function(e) NA)
    if (!is.na(cor_check) && abs(cor_check) < 0.999) {
      model_quad <- fit_model(~ mod_centered + mod_sq, dat_clean, use_random = n_studies >= 3)
      if (!is.null(model_quad)) {
        aic_quad <- tryCatch(AIC(model_quad), error = function(e) NA)
        p_quad   <- tryCatch(model_quad$pval[3], error = function(e) NA)
        quad_valid <- !is.na(p_quad) && p_quad < 0.05
        a <- tryCatch(model_quad$b[3], error = function(e) NA)
        b <- tryCatch(model_quad$b[2], error = function(e) NA)
        if (!is.na(a) && !is.na(b) && a != 0) {
          vertex <- -b / (2 * a) + x_mean
          vertex_type <- ifelse(a > 0, "minimum", "maximum")
        }
      }
    }
  }
  
  # Log
  model_log <- NULL; aic_log <- NA; p_log <- NA; log_valid <- FALSE
  if (moderator == "Total_Zn" && all(dat_clean$mod_transformed > 0, na.rm = TRUE)) {
    dat_clean$mod_log <- log(dat_clean$mod_transformed + 1)
    model_log <- fit_model(~ mod_log, dat_clean, use_random = n_studies >= 3)
    if (!is.null(model_log)) {
      aic_log <- tryCatch(AIC(model_log), error = function(e) NA)
      p_log   <- tryCatch(model_log$pval[2], error = function(e) NA)
      log_valid <- !is.na(p_log) && p_log < 0.05
    }
  }
  
  # Segmented
  model_seg <- NULL; aic_seg <- NA; p_seg <- NA; seg_valid <- FALSE; breakpoint <- NA
  if (nrow(dat_clean) >= 30) {
    bp_candidates <- tryCatch(quantile(dat_clean$mod_transformed, probs = seq(0.25, 0.75, 0.1)),
                              error = function(e) NULL)
    if (!is.null(bp_candidates) && length(bp_candidates) > 0) {
      best_model <- NULL; best_aic <- Inf; best_bp <- NA; best_p <- 1
      for (bp in bp_candidates) {
        dat_clean$seg_high <- as.numeric(dat_clean$mod_transformed > bp)
        dat_clean$mod_seg  <- dat_clean$mod_transformed * dat_clean$seg_high
        cor_check <- tryCatch(cor(dat_clean$mod_transformed, dat_clean$mod_seg,
                                  use = "complete.obs"), error = function(e) NA)
        if (!is.na(cor_check) && abs(cor_check) > 0.99) next
        temp <- fit_model(~ mod_transformed + mod_seg, dat_clean, use_random = n_studies >= 3)
        if (!is.null(temp) && length(temp$pval) >= 3) {
          p_int <- tryCatch(temp$pval[3], error = function(e) NA)
          if (!is.na(p_int) && p_int < 0.05) {
            aic_temp <- tryCatch(AIC(temp), error = function(e) Inf)
            if (!is.infinite(aic_temp) && aic_temp < best_aic) {
              best_aic <- aic_temp; best_model <- temp; best_bp <- bp; best_p <- p_int
            }
          }
        }
      }
      if (!is.null(best_model)) {
        model_seg <- best_model; aic_seg <- best_aic; p_seg <- best_p
        breakpoint <- best_bp; seg_valid <- TRUE
      }
    }
  }
  
  # Model selection
  models_list <- list()
  if (linear_valid && !is.null(model_linear))
    models_list[["Linear"]] <- list(model = model_linear, aic = aic_linear, p = p_linear, type = "Linear")
  if (quad_valid && !is.null(model_quad))
    models_list[["Quadratic"]] <- list(model = model_quad, aic = aic_quad, p = p_quad,
                                       type = "Quadratic", vertex = vertex, vertex_type = vertex_type)
  if (log_valid && !is.null(model_log))
    models_list[["Log"]] <- list(model = model_log, aic = aic_log, p = p_log, type = "Log")
  if (seg_valid && !is.null(model_seg))
    models_list[["Segmented"]] <- list(model = model_seg, aic = aic_seg, p = p_seg,
                                       type = "Segmented", breakpoint = breakpoint)
  
  if (length(models_list) == 0) {
    if (!is.null(model_linear)) {
      selected_model <- model_linear; model_type <- "Linear (non-significant)"; best_name <- "Linear_NS"
    } else return(NULL)
  } else {
    best_name <- names(models_list)[which.min(sapply(models_list, function(x) x$aic))]
    best <- models_list[[best_name]]
    selected_model <- best$model; model_type <- best$type
  }
  
  return(list(
    selected_model = selected_model, model_type = model_type, best_name = best_name,
    all_models = models_list, dat_clean = dat_clean,
    moderator = moderator, response_name = response_name, x_label = x_label,
    n = nrow(dat_clean), n_studies = n_studies,
    linear = list(p = p_linear, aic = aic_linear),
    quadratic = list(p = p_quad, aic = aic_quad, valid = quad_valid,
                     vertex = vertex, vertex_type = vertex_type),
    segmented = list(p = p_seg, aic = aic_seg, valid = seg_valid, breakpoint = breakpoint),
    log = list(p = p_log, aic = aic_log, valid = log_valid)
  ))
}

# ============================================
# 3. Re-run models to get results containing model objects
# ============================================
responses <- c(
  "Zn_shoot", "Zn_root", "Zn_plant",
  "Biomass_shoot", "Biomass_root", "Biomass_fruit",
  "Biomass_leave", "Plant_height",
  "BCF_shoot", "BCF_root", "BCF_plant",
  "P_shoot", "P_root"
)
moderators <- c("Total_Zn", "soil_pH")

all_results <- list()
for (resp in responses) {
  if (!resp %in% names(processed_data)) next
  dat <- processed_data[[resp]]
  if (is.null(dat)) next
  for (mod in moderators) {
    if (!mod %in% colnames(dat)) next
    cat(sprintf("🔄 Running model: %s ~ %s ... ", resp, mod))
    res <- compare_metareg_models(dat, mod, resp)
    if (!is.null(res)) {
      all_results[[paste(resp, mod, sep = "_")]] <- res
      cat("✅ (", res$model_type, ")\n", sep = "")
    } else cat("❌\n")
  }
}
cat("✅ Complete, total", length(all_results), "results\n")

# ============================================
# 4. Single panel plotting function
# ============================================
plot_one_panel <- function(res, tag = NULL) {
  dat   <- res$dat_clean
  mod   <- res$moderator
  x_var <- "mod_transformed"
  
  x_lab <- if (mod == "Total_Zn") expression("Total soil Zn (mg/kg, log"[10] * ")")
  else if (mod == "soil_pH") "Soil pH" else mod
  
  x_seq <- seq(min(dat[[x_var]], na.rm = TRUE),
               max(dat[[x_var]], na.rm = TRUE), length.out = 200)
  newdata <- data.frame(mod_transformed = x_seq)
  
  pred <- tryCatch({
    if (grepl("Quadratic", res$model_type)) {
      x_mean <- mean(dat$mod_transformed, na.rm = TRUE)
      newdata$mod_centered <- newdata$mod_transformed - x_mean
      newdata$mod_sq <- newdata$mod_centered^2
      predict(res$selected_model,
              newmods = cbind(newdata$mod_centered, newdata$mod_sq), addx = FALSE)
    } else if (grepl("Segmented", res$model_type)) {
      bp <- res$segmented$breakpoint
      newdata$seg_high <- as.numeric(newdata$mod_transformed > bp)
      newdata$mod_seg  <- newdata$mod_transformed * newdata$seg_high
      predict(res$selected_model,
              newmods = cbind(newdata$mod_transformed, newdata$mod_seg), addx = FALSE)
    } else if (grepl("Log", res$model_type)) {
      newdata$mod_log <- log(newdata$mod_transformed + 1)
      predict(res$selected_model, newmods = newdata$mod_log, addx = FALSE)
    } else {
      predict(res$selected_model, newmods = newdata$mod_transformed, addx = FALSE)
    }
  }, error = function(e) NULL)
  
  dat$prec <- 1 / dat$variance
  
  p <- ggplot() +
    geom_point(data = dat, aes(x = .data[[x_var]], y = lnRR, size = prec),
               alpha = 0.65, color = "#8B7355") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_size_continuous(
      range = c(1, 6),        # Point size range
      name = "Prec",          # Legend title
      limits = c(0, 300),     # Fixed range to ensure all panels unified (adjust to your max precision)
      breaks = c(50, 100, 150, 200) # Legend breaks
    ) +
    labs(x = x_lab, y = "lnRR (AMF effect)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold"))
  
  if (!is.null(pred)) {
    pred_df <- data.frame(x = x_seq, fit = pred$pred, se = pred$se)
    p <- p +
      geom_ribbon(data = pred_df,
                  aes(x = x, ymin = fit - 1.96*se, ymax = fit + 1.96*se),
                  fill = "#8B7355", alpha = 0.25) +
      geom_line(data = pred_df, aes(x = x, y = fit), color = "#5C4B2E", linewidth = 1)
  }
  
  if (grepl("Segmented", res$model_type) && !is.na(res$segmented$breakpoint))
    p <- p + geom_vline(xintercept = res$segmented$breakpoint,
                        linetype = "dotted", color = "steelblue")
  if (grepl("Quadratic", res$model_type) && !is.na(res$quadratic$vertex))
    p <- p + geom_vline(xintercept = res$quadratic$vertex,
                        linetype = "dotted", color = "steelblue")
  
  lab <- paste0("Model: ", res$model_type, "\n")
  if (grepl("Segmented", res$model_type) && !is.na(res$segmented$breakpoint)) {
    lab <- paste0(lab, "Breakpoint: ", round(res$segmented$breakpoint, 2), "\n",
                  "p = ", signif(res$segmented$p, 3), "\n")
  } else if (grepl("Quadratic", res$model_type) && !is.na(res$quadratic$vertex)) {
    lab <- paste0(lab, "p(quad) = ", signif(res$quadratic$p, 3), "\n",
                  "extreme at x = ", round(res$quadratic$vertex, 2), "\n")
  } else {
    lab <- paste0(lab, "p = ", signif(res$linear$p, 3), "\n")
  }
  lab <- paste0(lab, "n = ", res$n)
  
  p <- p + annotate("text", x = -Inf, y = -Inf, label = lab,
                    hjust = -0.05, vjust = -0.1, size = 3.2, color = "grey20")
  p <- p + ggtitle(paste(res$response_name, "~", mod))
  if (!is.null(tag)) p <- p + labs(tag = tag)
  p
}

# ============================================
# ============================================
# 5. Sort + combine pH and Soil Zn separately + output
# ============================================

# 1. Extract all results containing Total_Zn
zn_names <- names(all_results)[grepl("_Total_Zn$", names(all_results))]
# 2. Extract all results containing soil_pH
ph_names <- names(all_results)[grepl("_soil_pH$", names(all_results))]

# Predefined response variable priority order (to match reference figure layout)
preferred_order <- c(
  "Zn_shoot", "Zn_root", "Zn_plant", 
  "Biomass_shoot", "Biomass_root", "Biomass_fruit", "Biomass_leave", "Plant_height",
  "BCF_shoot", "BCF_root", "BCF_plant", "P_shoot", "P_root"
)

# Helper function: sort
sort_names <- function(names_vec, suffix) {
  # Remove suffix to get response variable names
  resp_names <- sub(paste0("_", suffix, "$"), "", names_vec)
  # Sort by priority order
  ordered_resp <- c(intersect(preferred_order, resp_names), 
                    setdiff(resp_names, preferred_order))
  return(paste0(ordered_resp, "_", suffix))
}

zn_names_ordered <- sort_names(zn_names, "Total_Zn")
ph_names_ordered <- sort_names(ph_names, "soil_pH")

# Define a function to plot and save
create_and_save_panel <- function(ordered_names, all_results, fig_name_prefix) {
  if (length(ordered_names) == 0) {
    cat("⚠️ No", fig_name_prefix, "results to plot, skipping.\n")
    return(NULL)
  }
  
  plot_list <- list()
  for (i in seq_along(ordered_names)) {
    tag <- if (i <= 26) paste0("(", letters[i], ")") else NULL
    plot_list[[i]] <- plot_one_panel(all_results[[ordered_names[i]]], tag)
  }
  
  # Compute rows and columns (max 4 columns)
  n_col <- min(4, length(plot_list))
  n_row <- ceiling(length(plot_list) / n_col)
  
  combined <- wrap_plots(plot_list, ncol = n_col) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "top",        # Legend at top
      legend.direction = "horizontal",# Horizontal arrangement
      legend.background = element_rect(fill = "grey95", color = NA), # Light gray background (reference figure style)
      legend.key = element_blank(),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )
  
  # Save PDF
  out_pdf <- file.path(FIGURES_DIR, paste0(fig_name_prefix, ".pdf"))
  ggsave(out_pdf, combined, width = 4.2 * n_col, height = 3.4 * n_row,
         dpi = 300, limitsize = FALSE)
  cat("✅ PDF saved:", out_pdf, "\n")
  
  # Save PNG
  out_png <- sub("\\.pdf$", ".png", out_pdf)
  ggsave(out_png, combined, width = 4.2 * n_col, height = 3.4 * n_row,
         dpi = 300, limitsize = FALSE)
  cat("✅ PNG saved:", out_png, "\n")
  
  return(combined)
}

# Plot and save separately
cat("\n--- Starting to plot Total Soil Zn panel ---\n")
fig_zn <- create_and_save_panel(zn_names_ordered, all_results, "Fig5A_MetaRegression_Total_Zn")

cat("\n--- Starting to plot Soil pH panel ---\n")
fig_ph <- create_and_save_panel(ph_names_ordered, all_results, "Fig5B_MetaRegression_Soil_pH")