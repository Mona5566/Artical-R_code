# ============================================================
# Low vs High Zn subgroup analysis forest plot (three-level random-effects model)
# Data source: all_processed_data.rds (list of 13 tibble)
# Fixes: 1) all pairwise comparison contrast matrix 2) FDR correction 3) I² statistics
#        4) UTF-8 encoding cleaning 5) three-level/two-level model traceability
#        6) smaller point shape 7) remove Zn group from y-axis 8) wider figure
# ============================================================

library(metafor)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- 0. Configuration ----

rds_path <- "D:/2026-03_oncemore/TryRun22-40/data/processed/all_processed_data.rds"

response_config <- list(
  Biomass_shoot = list(label = "Shoot Biomass"),
  Biomass_root  = list(label = "Root Biomass"),
  Biomass_fruit = list(label = "Fruit Biomass"),
  Biomass_leave = list(label = "Leaf Biomass"),
  Plant_height  = list(label = "Plant Height"),
  Zn_shoot      = list(label = "Zn Shoot"),
  Zn_root       = list(label = "Zn Root"),
  Zn_plant      = list(label = "Zn Plant"),
  BCF_shoot     = list(label = "BCF Shoot"),
  BCF_root      = list(label = "BCF Root"),
  BCF_plant     = list(label = "BCF Plant"),
  P_shoot       = list(label = "P Shoot"),
  P_root        = list(label = "P Root")
)

subgroup_vars   <- c("Type_soil","Type_AMF2","Duration","Type_root","Type_pH","Soil_P","Type_plant")
subgroup_labels <- c("Soil Type","AMF Type","Duration","Root Type","Soil pH","Soil P","Plant Type")

response_groups <- list(
  Biomass_Growth = list(
    vars   = c("Biomass_shoot","Biomass_root","Biomass_fruit","Biomass_leave","Plant_height"),
    labels = c("Shoot Biomass","Root Biomass","Fruit Biomass","Leaf Biomass","Plant Height"),
    ncols = 2, nrows = 3
  ),
  Zn_BCF = list(
    vars   = c("Zn_shoot","Zn_root","Zn_plant","BCF_shoot","BCF_root","BCF_plant"),
    labels = c("Zn Shoot","Zn Root","Zn Plant","BCF Shoot","BCF Root","BCF Plant"),
    ncols = 2, nrows = 3
  ),
  Phosphorus = list(
    vars   = c("P_shoot","P_root"),
    labels = c("P Shoot","P Root"),
    ncols = 2, nrows = 3
  )
)

output_dir <- "results/figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 0b. UTF-8 encoding cleaning utility ----

sanitize_utf8 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  bad <- !validUTF8(x)
  if (any(bad)) {
    for (enc in c("UTF-8","GBK","GB18030","latin1")) {
      y <- iconv(x, from = enc, to = "UTF-8", sub = "")
      if (!any(is.na(y[bad]))) { x <- y; break }
    }
    x <- iconv(x, from = "", to = "UTF-8", sub = "?")
  }
  x
}

# ---- 1. I² statistic (report model level) ----

compute_I2 <- function(fit, dat_clean) {
  sigma2_vals <- fit$sigma2
  n_varcomp   <- length(sigma2_vals)
  
  if (n_varcomp < 2) {
    sigma2_total <- sum(sigma2_vals)
    v_typical <- median(dat_clean$variance, na.rm = TRUE)
    return(data.frame(
      I2_total   = sigma2_total / (sigma2_total + v_typical) * 100,
      I2_between = NA_real_,
      I2_within  = NA_real_,
      n_varcomp  = n_varcomp,
      model_type = "2-level (fallback)"
    ))
  }
  
  s_b <- sigma2_vals[1]; s_w <- sigma2_vals[2]
  v_typical <- median(dat_clean$variance, na.rm = TRUE)
  tot <- s_b + s_w + v_typical
  data.frame(
    I2_total   = (s_b + s_w) / tot * 100,
    I2_between = s_b / tot * 100,
    I2_within  = s_w / tot * 100,
    n_varcomp  = n_varcomp,
    model_type = "3-level"
  )
}

# ---- 2. All pairwise comparison contrast matrix ----

make_pairwise_contrasts <- function(p) {
  pairs <- combn(p, 2); n_pairs <- ncol(pairs)
  L_mat <- matrix(0, nrow = n_pairs, ncol = p)
  for (k in 1:n_pairs) { L_mat[k, pairs[1,k]] <- 1; L_mat[k, pairs[2,k]] <- -1 }
  L_mat
}

# ---- 3. Three-level random-effects model subgroup analysis ----

run_subgroup_analysis <- function(dat, var, cat_label, min_n = 3) {
  if (!var %in% colnames(dat)) return(NULL)
  
  dat[[var]] <- sanitize_utf8(dat[[var]])
  dat_clean <- dat[complete.cases(dat[[var]], dat$lnRR, dat$variance), ]
  dat_clean <- dat_clean[is.finite(dat_clean$lnRR) & dat_clean$variance > 0, ]
  if (nrow(dat_clean) < min_n) return(NULL)
  
  dat_clean[[var]] <- trimws(as.character(dat_clean[[var]]))
  dat_clean <- dat_clean[dat_clean[[var]] != "" & dat_clean[[var]] != "NA", ]
  groups_levels <- sort(unique(dat_clean[[var]]))
  if (length(groups_levels) < 2) return(NULL)
  counts <- table(dat_clean[[var]])
  if (any(counts < min_n)) return(NULL)
  
  if (!"effect_id" %in% colnames(dat_clean)) {
    dat_clean <- dat_clean %>% group_by(study_id) %>% mutate(effect_id = row_number()) %>% ungroup()
  }
  
  mod_formula <- as.formula(paste0("~ 0 + ", var))
  fit <- tryCatch({
    rma.mv(yi = lnRR, V = variance, data = dat_clean,
           random = ~ 1 | study_id/effect_id, mods = mod_formula, method = "REML",
           control = list(optimizer = "optim", optmethod = "Nelder-Mead", rel.tol = 1e-8))
  }, error = function(e) {
    tryCatch({
      rma.mv(yi = lnRR, V = variance, data = dat_clean,
             random = ~ 1 | study_id, mods = mod_formula, method = "REML",
             control = list(optimizer = "optim", optmethod = "Nelder-Mead"))
    }, error = function(e2) NULL)
  })
  if (is.null(fit)) return(NULL)
  
  coef_sum <- coef(summary(fit))
  rn <- rownames(coef_sum)
  group_names <- ifelse(startsWith(rn, var),
                        substring(rn, nchar(var) + 1L),
                        rn)
  group_names <- trimws(group_names)
  
  n_vals <- as.numeric(counts[group_names]); n_vals[is.na(n_vals)] <- min(counts)
  
  QM_val <- NA; QMp_val <- NA; pairwise_results <- NULL
  if (length(groups_levels) > 1) {
    L_mat <- make_pairwise_contrasts(length(groups_levels)); colnames(L_mat) <- groups_levels
    qm <- tryCatch({ t <- anova(fit, L = L_mat); list(QM = t$QM, QMp = t$QMp) },
                   error = function(e) list(QM = NA, QMp = NA))
    QM_val <- qm$QM; QMp_val <- qm$QMp
    pairwise_list <- list(); pairs <- combn(length(groups_levels), 2)
    for (k in 1:ncol(pairs)) {
      L_row <- matrix(L_mat[k, ], nrow = 1)
      pw <- tryCatch({
        res <- anova(fit, L = L_row)
        data.frame(group1 = groups_levels[pairs[1,k]], group2 = groups_levels[pairs[2,k]],
                   diff = res$QM, pval = res$QMp, stringsAsFactors = FALSE)
      }, error = function(e) NULL)
      if (!is.null(pw)) pairwise_list[[length(pairwise_list)+1]] <- pw
    }
    if (length(pairwise_list) > 0) {
      pairwise_results <- bind_rows(pairwise_list)
      pairwise_results$pval_fdr <- p.adjust(pairwise_results$pval, method = "BH")
    }
  }
  
  I2_df <- compute_I2(fit, dat_clean)
  df <- data.frame(
    category   = cat_label, group = group_names,
    rr         = coef_sum[, "estimate"],
    lower      = coef_sum[, "ci.lb"],
    upper      = coef_sum[, "ci.ub"],
    n          = n_vals,
    p_value    = coef_sum[, "pval"],
    tau2       = fit$tau2,
    I2_total   = I2_df$I2_total,
    I2_between = I2_df$I2_between,
    I2_within  = I2_df$I2_within,
    n_varcomp  = I2_df$n_varcomp,
    model_type = I2_df$model_type,
    QM         = QM_val,
    QMp        = QMp_val,
    stringsAsFactors = FALSE
  )
  df$p_lab <- ifelse(df$p_value < 0.001, "P<0.001",
                     ifelse(df$p_value < 0.01,  "P<0.01",
                            ifelse(df$p_value < 0.05,  "P<0.05",
                                   paste0("P=", round(df$p_value, 3)))))
  df$signif <- ifelse(df$p_value < 0.001, "***",
                      ifelse(df$p_value < 0.01,  "**",
                             ifelse(df$p_value < 0.05,  "*", "")))
  attr(df, "pairwise") <- pairwise_results
  df
}

# ---- 4. Read RDS (list of 13 tibble), clean all columns UTF-8 ----

cat("Reading RDS data...\n")
if (!file.exists(rds_path)) stop("RDS file not found: ", rds_path)
rds_obj <- readRDS(rds_path)
stopifnot(is.list(rds_obj))
cat("RDS list length: ", length(rds_obj), "\n")
cat("Names: ", paste(names(rds_obj), collapse = ", "), "\n")

for (resp_name in names(rds_obj)) {
  df0 <- rds_obj[[resp_name]]
  if (!is.data.frame(df0)) next
  for (col in colnames(df0)) {
    if (is.character(df0[[col]]) || is.factor(df0[[col]])) {
      x <- as.character(df0[[col]])
      bad <- !validUTF8(x)
      if (any(bad)) {
        cat(sprintf("  ⚠ Non-UTF-8: %s$%s, first bad value: %s\n",
                    resp_name, col, x[which(bad)[1]]))
      }
    }
  }
}

rds_datasets <- list()
for (resp_name in names(response_config)) {
  df <- rds_obj[[resp_name]]
  if (is.null(df)) { cat(sprintf("  ⚠ RDS missing element: %s\n", resp_name)); next }
  
  for (col in colnames(df)) {
    if (is.character(df[[col]]) || is.factor(df[[col]])) {
      df[[col]] <- sanitize_utf8(df[[col]])
    }
  }
  
  if ("se" %in% colnames(df) && !("variance" %in% colnames(df))) df$variance <- df$se^2
  df$Zn_Level <- tolower(trimws(as.character(df$Zn_Level)))
  
  need <- c("lnRR","variance","Zn_Level","study_id")
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) {
    cat(sprintf("  ⚠ %s missing columns: %s, skipping\n", resp_name, paste(miss, collapse=", ")))
    next
  }
  rds_datasets[[resp_name]] <- df
}

# ---- 5. Run models ----

cat("\nStarting three-level random-effects model subgroup analysis...\n")
all_results <- list(); all_pairwise <- list()

for (resp_name in names(response_config)) {
  cfg <- response_config[[resp_name]]
  dat <- rds_datasets[[resp_name]]
  if (is.null(dat)) { cat(sprintf("  ⚠ Missing data: %s\n", resp_name)); next }
  if (nrow(dat) < 3)  { cat(sprintf("  ⚠ Insufficient rows: %s\n", resp_name)); next }
  
  all_results[[resp_name]] <- list()
  
  for (zn_level in c("low","high")) {
    dat_sub <- dat[dat$Zn_Level == zn_level, ]
    if (nrow(dat_sub) < 3) next
    
    result_list <- list()
    for (j in seq_along(subgroup_vars)) {
      result <- run_subgroup_analysis(dat_sub, subgroup_vars[j], subgroup_labels[j])
      if (!is.null(result)) {
        pw <- attr(result, "pairwise")
        if (!is.null(pw)) {
          pw$response <- cfg$label; pw$resp_var <- resp_name
          pw$zn_group <- ifelse(zn_level == "low", "Low Zn", "High Zn")
          all_pairwise[[length(all_pairwise)+1]] <- pw
        }
        result_list[[j]] <- result
      }
    }
    
    if (length(result_list) > 0) {
      combined <- bind_rows(result_list)
      combined$zn_group <- ifelse(zn_level == "low", "Low Zn", "High Zn")
      combined$response <- cfg$label
      combined$resp_var <- resp_name
      all_results[[resp_name]][[zn_level]] <- combined
      cat(sprintf("  ✓ %s (%s): %d subgroups\n", cfg$label, zn_level, length(result_list)))
      
      n_fallback <- sum(combined$model_type == "2-level (fallback)", na.rm = TRUE)
      if (n_fallback > 0) {
        fb_cats <- unique(combined$category[combined$model_type == "2-level (fallback)"])
        cat(sprintf("      ⚠ Of which %d/%d subgroups degenerated to two-level model: %s\n",
                    n_fallback, nrow(combined), paste(fb_cats, collapse = ", ")))
      }
    }
  }
}

# ---- 6. Compact forest plot for single response variable ----
make_compact_forest <- function(resp_name, resp_label, all_results,
                                global_xmin, global_xmax) {
  low_df  <- all_results[[resp_name]][["low"]]
  high_df <- all_results[[resp_name]][["high"]]
  if (is.null(low_df) && is.null(high_df)) return(NULL)
  
  rows <- list()
  if (!is.null(low_df))  rows[[length(rows) + 1]] <- low_df
  if (!is.null(high_df)) rows[[length(rows) + 1]] <- high_df
  combined <- bind_rows(rows)
  
  combined$category <- factor(combined$category, levels = subgroup_labels)
  
  # ★ Forced sorting: by category, subgroup name, Low before High
  combined <- combined %>%
    mutate(category_num = as.numeric(category)) %>%
    arrange(category_num, group, desc(zn_group == "Low Zn"))
  
  # ---- Y-axis mapping ----
  combined$y_key <- paste0(combined$group, "|", combined$zn_group)
  y_levels <- unique(combined$y_key)
  combined$y_label <- factor(combined$y_key, levels = rev(y_levels))
  combined$y_num <- as.numeric(combined$y_label)
  
  axis_map <- combined %>%
    select(y_num, group) %>%
    distinct() %>%
    arrange(y_num)
  
  y_breaks <- axis_map$y_num
  y_labels <- axis_map$group
  
  # Compute canvas base bounds and text positions
  x_span <- global_xmax - global_xmin
  if(x_span == 0) x_span <- 1
  
  x_left_text  <- global_xmin - x_span * 0.05
  x_right_text <- global_xmax + x_span * 0.05
  
  # ★ Generate independent text data frame (text_df) from main dataset
  combined <- combined %>%
    mutate(
      p_lab2 = ifelse(is.na(p_value), "P=NA",
                      ifelse(p_value < 0.001, "P<0.001",
                             ifelse(p_value < 0.01,  "P<0.01",
                                    ifelse(p_value < 0.05,  "P<0.05",
                                           paste0("P=", sprintf("%.3f", p_value)))))), 
      n_pad = sprintf("%-5s", paste0("n=", n)), 
      label_txt = paste0(n_pad, " ", p_lab2),
      signif_txt = ifelse(is.na(signif), "", signif)
    )
  
  # Text layer data, fixed x, y and label to prevent misalignment
  text_df <- combined %>%
    mutate(
      x_pos = ifelse(zn_group == "Low Zn", x_left_text, x_right_text),
      hjust_v = ifelse(zn_group == "Low Zn", 1, 0)
    )
  
  ggplot(combined, aes(x = rr, y = y_num, color = zn_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15, linewidth = 0.5) +
    geom_point(aes(shape = zn_group, fill = zn_group), size = 1.2) +
    
    # ★ Stars (bound directly to combined, positioned by upper)
    geom_text(data = text_df, 
              aes(x = upper + x_span * 0.015, y = y_num, label = signif_txt, color = zn_group), 
              hjust = 0, vjust = 0.5, size = 3.8, fontface = "bold",
              show.legend = FALSE) +
    
    # ★ Key fix: both sides text use independent dataset and explicitly pass hjust
    geom_text(data = text_df,
              aes(x = x_pos, y = y_num, label = label_txt, color = zn_group),
              hjust = text_df$hjust_v, vjust = 0.5, size = 3.0, family = "mono",
              show.legend = FALSE) +
    
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    scale_color_manual(values = c("Low Zn" = "#0072B2", "High Zn" = "#D55E00"), name = "Zn Level") +
    scale_fill_manual(values  = c("Low Zn" = "#0072B2", "High Zn" = "#D55E00"), name = "Zn Level") +
    scale_shape_manual(values = c("Low Zn" = 22, "High Zn" = 24), name = "Zn Level") +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.85, 0.85))
    ) +
    
    scale_y_continuous(
      breaks = y_breaks, 
      labels = y_labels,
      expand = expansion(add = 1.0)
    ) +
    
    coord_cartesian(clip = "off") +
    
    labs(x = "lnRR", y = NULL, title = resp_label) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#555555", linewidth = 0.5),
      strip.background = element_rect(fill = "#E6EDF2", color = NA),
      
      strip.text.y = element_text(size = 9.5, face = "bold.italic", color = "#2C5F8A", angle = 0),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "#1A1A1A"),
      
      legend.position = "none",
      plot.margin = margin(top = 10, right = 45, bottom = 10, left = 10, unit = "pt")
    )
}
# ---- 7. Generate A4 portrait grouped PDF ----

cat("\nGenerating A4 portrait grouped PDF...\n")
for (group_name in names(response_groups)) {
  gc <- response_groups[[group_name]]
  all_vals <- c()
  for (rn in gc$vars) for (zl in c("low","high")) {
    df <- all_results[[rn]][[zl]]
    if (!is.null(df)) all_vals <- c(all_vals, df$lower, df$upper)
  }
  all_vals <- all_vals[is.finite(all_vals)]
  if (length(all_vals) == 0) { cat(sprintf("  ⚠ %s: no available data\n", group_name)); next }
  global_xmin <- min(min(all_vals) - 0.3, -0.5)
  global_xmax <- max(max(all_vals) + 0.3,  0.5)
  
  plot_list <- list()
  for (i in seq_along(gc$vars)) {
    p <- make_compact_forest(gc$vars[i], gc$labels[i], all_results, global_xmin, global_xmax)
    if (!is.null(p)) plot_list[[i]] <- p
  }
  if (length(plot_list) == 0) { cat(sprintf("  ⚠ %s: no available plots\n", group_name)); next }
  
  # ★ Key modification: add left/right spacing when combining
  p_final <- wrap_plots(plot_list, ncol = gc$ncols, nrow = gc$nrows) +
    plot_annotation(
      title = paste0("Low Zn vs High Zn: ", group_name),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "#1A1A1A"))
    )
  
  # ★ Widen A4 portrait to 10 x 14 to give more space to plots
  pdf_path <- file.path(output_dir, paste0("low_high_subgroup_", group_name, ".pdf"))
  ggsave(pdf_path, p_final, width = 10, height = 16, device = "pdf", limitsize = FALSE)
  cat(sprintf("  ✓ %s: %s\n", group_name, pdf_path))
}
# ---- 8. Save summary data ----

all_combined <- list()
for (rn in names(all_results)) for (zl in names(all_results[[rn]]))
  all_combined[[length(all_combined)+1]] <- all_results[[rn]][[zl]]
if (length(all_combined) > 0) {
  all_df <- bind_rows(all_combined)
  write.csv(all_df, file.path(output_dir, "all_subgroup_results.csv"), row.names = FALSE)
  cat(sprintf("✓ Summary data: %d rows\n", nrow(all_df)))
  
  cat("\n---- Model level summary ----\n")
  print(table(all_df$model_type, useNA = "ifany"))
  if (any(all_df$model_type == "2-level (fallback)", na.rm = TRUE)) {
    fb <- all_df %>%
      filter(model_type == "2-level (fallback)") %>%
      distinct(resp_var, zn_group, category) %>%
      arrange(resp_var, zn_group, category)
    cat("\nDegenerated to two-level model (response variable, Zn level, subgroup):\n")
    print(fb, n = nrow(fb))
  }
}

if (length(all_pairwise) > 0) {
  pw_df <- bind_rows(all_pairwise)
  write.csv(pw_df, file.path(output_dir, "pairwise_comparisons.csv"), row.names = FALSE)
  cat(sprintf("\n✓ Pairwise comparison results: %d pairs, with FDR correction\n", nrow(pw_df)))
  cat(sprintf("  Significant proportion before FDR correction: %.1f%%\n",
              mean(pw_df$pval     < 0.05, na.rm = TRUE) * 100))
  cat(sprintf("  Significant proportion after FDR correction: %.1f%%\n",
              mean(pw_df$pval_fdr < 0.05, na.rm = TRUE) * 100))
}

cat("\n✅ Analysis complete!\n")