# ============================================
# Module 12: Grouped variable importance analysis (supports Zn_Level grouping)
# Scheme A: Ecology meta random forest standard paradigm
# Rule: normalize within single model; average across response variables, no second normalization, allow sum != 100%
# Support: All / Low Zn / High Zn three groups analyzed separately
# [Fixed version] Fixed theme duplicate parameter error, color bug, chart text ambiguity, added diagnostic checks
# ============================================

source("00_config.R")

library(randomForest)
library(ggplot2)
library(dplyr)
library(forcats)
library(patchwork)
library(tidyr)
library(readr)

# ============================================
# 1. Create output directory
# ============================================
GROUPED_IMP_DIR <- file.path(FIGURES_DIR, "grouped_importance")
dir.create(GROUPED_IMP_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================
# 2. Color mapping
# ============================================
GROUP_COLORS <- c(
  "Soil properties" = "#f2a65a",
  "AMF traits" = "#5d87c3",
  "Environment" = "#88bc64",
  "Plant traits" = "#bc7c9c"
)

fill_colors <- GROUP_COLORS
fill_colors[["Other"]] <- "gray70"
GROUP_COLORS_ALL <- c(GROUP_COLORS, "Other" = "gray70") # Full version with Other

# Zn_Level colors
ZN_LEVEL_COLORS <- c("low" = "#2E86AB", "high" = "#A23B72")

# ============================================
# 3. Variable grouping
# ============================================
VAR_GROUPS <- list(
  "Soil properties" = c("soil_pH", "Total_Zn", "Soil_P", "Type_soil"),
  "AMF traits" = c("Type_AMF2"),
  "Environment" = c("Duration"),
  "Plant traits" = c("Type_root", "Type_plant")
)

# ============================================
# 4. Core: run grouped random forest (supports group labels)
# ============================================
run_grouped_rf <- function(dat, response_name, zn_level = NULL, group_label = "All") {
  
  cat(sprintf("\n=== %s (%s) ===\n", response_name, group_label))
  
  available_vars <- c()
  for(group in names(VAR_GROUPS)) {
    for(var in VAR_GROUPS[[group]]) {
      if(var %in% colnames(dat)) {
        available_vars <- c(available_vars, var)
      }
    }
  }
  
  if(length(available_vars) < 2) {
    cat(sprintf("  Insufficient available variables (%d), skipping\n", length(available_vars)))
    return(NULL)
  }
  
  cat(sprintf("  Available variables: %s\n", paste(available_vars, collapse = ", ")))
  
  rf_data <- dat[, c("lnRR", available_vars), drop = FALSE]
  
  for(var in available_vars) {
    if(is.character(rf_data[[var]]) || is.factor(rf_data[[var]])) {
      rf_data[[var]] <- as.factor(rf_data[[var]])
    } else if(is.numeric(rf_data[[var]])) {
      if(var %in% c("Soil_P", "Duration", "Type_soil", "Type_AMF2", "Type_root", "Type_plant")) {
        rf_data[[var]] <- as.factor(rf_data[[var]])
      }
    }
  }
  
  rf_data <- na.omit(rf_data)
  
  if(nrow(rf_data) < 20) {
    cat(sprintf("  Insufficient sample size (n=%d), skipping\n", nrow(rf_data)))
    return(NULL)
  }
  
  cat(sprintf("  Final sample: n=%d, features=%d\n", nrow(rf_data), ncol(rf_data)-1))
  
  set.seed(123)
  rf_model <- randomForest(
    lnRR ~ .,
    data = rf_data,
    importance = TRUE,
    ntree = 500,
    mtry = max(1, floor(sqrt(ncol(rf_data) - 1)))
  )
  
  r2 <- tail(rf_model$rsq, 1)
  if(is.null(r2) || length(r2) == 0 || is.na(r2)) {
    r2 <- 0
  }
  explain_R2 <- round(r2 * 100, 1)
  
  imp_matrix <- rf_model$importance
  
  if("%IncMSE" %in% colnames(imp_matrix)) {
    imp_raw <- imp_matrix[, "%IncMSE"]
  } else if("IncNodePurity" %in% colnames(imp_matrix)) {
    imp_raw <- imp_matrix[, "IncNodePurity"]
  } else {
    if(ncol(imp_matrix) == 1) {
      imp_raw <- imp_matrix[, 1]
    } else {
      cat("  Cannot identify importance column\n")
      return(NULL)
    }
  }
  
  imp_df <- data.frame(
    Variable = names(imp_raw),
    Importance_raw = as.numeric(imp_raw)
  ) %>% arrange(desc(Importance_raw))
  
  cat(sprintf("  Model R² = %.1f%%\n", explain_R2))
  cat("  Importance ranking:\n")
  print(head(imp_df, 10))
  
  imp_df <- imp_df %>%
    mutate(
      Group = case_when(
        Variable %in% VAR_GROUPS[["Soil properties"]] ~ "Soil properties",
        Variable %in% VAR_GROUPS[["AMF traits"]] ~ "AMF traits",
        Variable %in% VAR_GROUPS[["Environment"]] ~ "Environment",
        Variable %in% VAR_GROUPS[["Plant traits"]] ~ "Plant traits",
        TRUE ~ "Other"
      )
    ) %>%
    filter(Importance_raw > 0) # Remove negative importance variables (no predictive value)
  
  if(nrow(imp_df) == 0) {
    cat("  No valid importance values\n")
    return(NULL)
  }
  
  # Normalize within single model: sum ≈ 100% for same response variable
  imp_df <- imp_df %>%
    mutate(
      Importance_pct = Importance_raw / sum(Importance_raw) * 100
    )
  
  cat("  Group contribution (normalized %, sum ≈ 100% within model; Other not shown in console):\n")
  group_summary <- imp_df %>%
    filter(Group != "Other") %>%
    group_by(Group) %>%
    summarise(Total = sum(Importance_pct), .groups = "drop") %>%
    arrange(desc(Total))
  print(group_summary)
  
  return(list(
    response = response_name,
    zn_level = zn_level,
    group_label = group_label,
    importance = imp_df,
    n_obs = nrow(rf_data),
    r_squared = explain_R2,
    n_vars = ncol(rf_data) - 1,
    model = rf_model
  ))
}

# ============================================
# 5. Plotting function (supports grouping)
# ============================================
plot_grouped_importance <- function(result, zn_level = NULL) {
  
  imp_df <- result$importance
  
  # Group colors
  if(!is.null(zn_level) && zn_level %in% names(ZN_LEVEL_COLORS)) {
    main_color <- ZN_LEVEL_COLORS[zn_level]
  } else {
    main_color <- "#4C72B0"
  }
  
  # Left panel: donut chart
  group_data <- imp_df %>%
    filter(Group != "Other") %>%
    group_by(Group) %>%
    summarise(Total = sum(Importance_pct), .groups = "drop") %>%
    mutate(
      Group = factor(Group, levels = names(GROUP_COLORS)),
      label = paste0(round(Total, 1), "%"),
      ypos = cumsum(Total) - 0.5 * Total
    )
  
  has_groups <- nrow(group_data) > 0
  
  if(has_groups) {
    p1 <- ggplot(group_data, aes(x = 2, y = Total, fill = Group)) +
      geom_bar(stat = "identity", width = 1, alpha = 0.9, color = "white", linewidth = 0.5) +
      coord_polar(theta = "y", start = 0) +
      xlim(0.5, 2.5) +
      geom_text(aes(label = label, y = ypos), 
                size = 4.5, fontface = "bold", color = "white") +
      annotate("text", x = 0, y = 0, 
               label = "Within-model\ncontribution",
               size = 4.5, fontface = "bold", color = "gray30") +
      scale_fill_manual(values = GROUP_COLORS) +
      labs(
        title = paste("Relative Contribution by Group -", result$group_label),
        subtitle = "Within-model normalized %IncMSE (sum to 100% for this response)"
      ) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 8, color = "gray40"),
        legend.position = "none"
      )
  } else {
    p1 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No grouped variables", size = 5) +
      theme_void()
  }
  
  # Right panel: bar chart
  imp_df <- imp_df %>%
    mutate(Variable = factor(Variable, levels = rev(Variable)))
  
  p2 <- ggplot(imp_df, aes(x = Importance_pct, y = Variable, fill = Group)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_text(aes(label = round(Importance_pct, 1)), 
              hjust = -0.15, size = 3.2) +
    scale_fill_manual(values = fill_colors, name = "Variable Group") +
    labs(
      x = "Relative importance of predictors (%)",
      y = NULL
    ) +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("Explained variation: ", result$r_squared, "%"),
             hjust = 1.05, vjust = -0.5, size = 3.8, 
             fontface = "italic", color = "gray30") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      plot.margin = unit(c(0.5, 0.3, 0.3, 0.3), "cm")
    )
  
  # Combine plots
  combined <- p1 + p2 + 
    plot_layout(widths = c(0.3, 0.7)) +
    plot_annotation(
      title = paste(result$response, "-", result$group_label),
      subtitle = paste0("n = ", result$n_obs, " observations, ", 
                        result$n_vars, " predictors"),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40")
      )
    )
  
  return(combined)
}

# ============================================
# 6. Batch run (supports grouping)
# ============================================
run_all_grouped_importance <- function(processed_data, zn_level = NULL, group_label = "All") {
  
  cat("\n")
  cat("========================================\n")
  cat(sprintf("  Grouped variable importance analysis - %s\n", group_label))
  cat("  Importance metric: %IncMSE | Scheme A: normalize within model, cross-response mean not rescaled\n")
  cat("  Style: donut chart + bar chart\n")
  cat("========================================\n")
  cat(sprintf("Number of response variables: %d\n", length(processed_data)))
  
  # Create grouped output directory
  if(is.null(zn_level)) {
    out_dir <- file.path(GROUPED_IMP_DIR, "all")
  } else {
    out_dir <- file.path(GROUPED_IMP_DIR, zn_level)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("Output directory: %s\n", out_dir))
  cat("\n")
  
  all_results <- list()
  
  for(resp in names(processed_data)) {
    dat <- processed_data[[resp]]
    if(is.null(dat)) next
    
    result <- tryCatch({
      run_grouped_rf(dat, resp, zn_level, group_label)
    }, error = function(e) {
      cat(sprintf("  Error: %s\n", e$message))
      return(NULL)
    })
    
    if(!is.null(result)) {
      all_results[[resp]] <- result
      
      p <- plot_grouped_importance(result, zn_level)
      
      safe_name <- gsub("[^A-Za-z0-9_]", "_", resp)
      
      # Save to grouped directory
      ggsave(file.path(out_dir, paste0(safe_name, "_grouped_importance.png")),
             p, width = 11, height = 6.5, dpi = 350, bg = "white")
      ggsave(file.path(out_dir, paste0(safe_name, "_grouped_importance.pdf")),
             p, width = 11, height = 6.5, device = "pdf")
      
      cat(sprintf("  ✓ %s: R²=%.1f%%, variables=%d, n=%d\n", 
                  resp, result$r_squared, result$n_vars, result$n_obs))
    }
  }
  
  # ---- Summary plot (Scheme A: remove second normalization!!) ----
  if(length(all_results) > 1) {
    cat("\n--- Generating summary plot (all response variables) ---\n")
    
    # 1. Merge all importance data
    all_imp <- bind_rows(lapply(names(all_results), function(resp) {
      df <- all_results[[resp]]$importance
      df$Response <- resp
      return(df[, c("Variable", "Group", "Importance_pct", "Response")])
    }))
    
    cat(sprintf("Merged data: %d rows\n", nrow(all_imp)))
    cat(sprintf("Unique variables: %d\n", length(unique(all_imp$Variable))))
    cat(sprintf("Unique response variables: %d\n", length(unique(all_imp$Response))))
    
    # 2. Check variable grouping consistency
    cat("\n--- Checking variable grouping consistency ---\n")
    
    var_group_check <- all_imp %>%
      group_by(Variable) %>%
      summarise(
        groups = paste(unique(Group), collapse = ", "),
        n_groups = n_distinct(Group),
        .groups = "drop"
      ) %>%
      filter(n_groups > 1)
    
    if(nrow(var_group_check) > 0) {
      cat("Found inconsistent variable grouping:\n")
      print(var_group_check)
      
      # Majority vote to unify variable grouping
      cat("\nUsing majority vote to fix grouping...\n")
      variable_group_mapping <- all_imp %>%
        group_by(Variable, Group) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(Variable) %>%
        slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
        ungroup()
      
      variable_group_mapping <- variable_group_mapping[, c("Variable", "Group")]
      colnames(variable_group_mapping) <- c("Variable", "Consistent_Group")
      
      all_imp_clean <- all_imp %>%
        dplyr::select(-Group) %>%
        left_join(variable_group_mapping, by = "Variable") %>%
        rename(Group = Consistent_Group)
      
      cat("Grouping fix complete\n")
    } else {
      cat("All variable groupings consistent ✓\n")
      all_imp_clean <- all_imp
    }
    
    # 3. Fill missing variables with 0 (no positive importance for that response variable)
    all_vars <- unique(all_imp_clean$Variable)
    all_responses <- names(all_results)
    
    cat(sprintf("\nCreating complete data matrix: %d variables × %d response variables\n", 
                length(all_vars), length(all_responses)))
    
    complete_imp <- expand.grid(
      Variable = all_vars,
      Response = all_responses,
      stringsAsFactors = FALSE
    ) %>%
      left_join(all_imp_clean, by = c("Variable", "Response")) %>%
      mutate(
        Importance_pct = ifelse(is.na(Importance_pct), 0, Importance_pct),
        Group = ifelse(is.na(Group), "Other", Group)
      )
    
    # 4. Compute average importance
    avg_imp <- complete_imp %>%
      group_by(Variable, Group) %>%
      summarise(
        Mean_Pct = mean(Importance_pct, na.rm = TRUE),
        n_responses = sum(Importance_pct > 0),
        .groups = "drop"
      ) %>%
      arrange(desc(Mean_Pct))
    
    # ======================
    # New: diagnostic check code
    # ======================
    cat("\n==== 【Diagnostic Check】Sum of Mean_Pct for all variables ====\n")
    total_all_vars <- sum(avg_imp$Mean_Pct)
    cat(sprintf("Sum of all variable means = %.3f \n", total_all_vars))
    if(abs(total_all_vars - 100) > 0.15){
      warning("⚠️ Note: cross-response mean sum need not equal 100; this paradigm normalizes within single model, cross-response mean is not re-standardized")
    }
    # ======================
    
    cat(sprintf("avg_imp: %d rows, %d columns\n", nrow(avg_imp), ncol(avg_imp)))
    
    # Save summary table
    suffix <- ifelse(is.null(zn_level), "all", zn_level)
    write_csv(avg_imp, file.path(TABLES_DIR, paste0("grouped_importance_summary_", suffix, ".csv")))
    
    # 5. Prepare donut chart data (including Other group)
    avg_group <- avg_imp %>%
      group_by(Group) %>%
      summarise(Total = sum(Mean_Pct), .groups = "drop") %>%
      mutate(
        Group = factor(Group, levels = c(names(GROUP_COLORS), "Other")),
        label = paste0(round(Total, 1), "%"),
        ypos = cumsum(Total) - 0.5 * Total
      )
    
    # 6. Draw donut chart [Fix: use GROUP_COLORS_ALL]
    if(nrow(avg_group) > 0) {
      p_avg_donut <- ggplot(avg_group, aes(x = 2, y = Total, fill = Group)) +
        geom_bar(stat = "identity", width = 1, alpha = 0.9, color = "white", linewidth = 0.5) +
        coord_polar(theta = "y", start = 0) +
        xlim(0.5, 2.5) +
        geom_text(aes(label = label, y = ypos), 
                  size = 4.5, fontface = "bold", color = "white") +
        annotate("text", x = 0, y = 0, 
                 label = "Mean importance\nacross responses",
                 size = 4.5, fontface = "bold", color = "gray30") +
        scale_fill_manual(values = GROUP_COLORS_ALL) + # Fix color bug
        labs(
          title = paste("Average Contribution by Group -", group_label),
          subtitle = paste0("Mean normalized %IncMSE across predictors; values do NOT enforce sum to 100%")
        ) +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
          legend.position = "none"
        )
    } else {
      p_avg_donut <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No grouped variables", size = 5) +
        theme_void()
    }
    
    # 7. Draw bar chart [removed duplicate panel.grid.minor parameter!]
    avg_imp_for_bar <- avg_imp %>%
      arrange(desc(Mean_Pct))
    
    unique_vars <- unique(avg_imp_for_bar$Variable)
    
    p_avg_bar <- avg_imp_for_bar %>%
      mutate(
        Variable = factor(Variable, levels = unique_vars)
      ) %>%
      ggplot(aes(x = Mean_Pct, y = Variable, fill = Group)) +
      geom_col(width = 0.7, alpha = 0.85) +
      geom_text(aes(label = round(Mean_Pct, 1)), 
                hjust = -0.15, size = 3) +
      scale_fill_manual(values = fill_colors, name = "Variable Group") +
      labs(
        x = "Mean relative importance of predictors (%)",
        y = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9)
      )
    
    # 8. Combine summary plot
    p_avg_combined <- p_avg_donut + p_avg_bar + 
      plot_layout(widths = c(0.3, 0.7)) +
      plot_annotation(
        title = paste("Overall Variable Importance -", group_label),
        subtitle = paste0("Mean normalized %IncMSE across ", length(all_results),
                          " response variables"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40")
        )
      )
    
    # 9. Save
    ggsave(file.path(out_dir, "Average_Grouped_Importance.png"),
           p_avg_combined, width = 12, height = 8, dpi = 350, bg = "white")
    ggsave(file.path(out_dir, "Average_Grouped_Importance.pdf"),
           p_avg_combined, width = 12, height = 8, device = "pdf")
    
    cat("  ✓ Summary plot saved\n")
    
    # 10. Output average contribution by group
    group_summary <- avg_imp %>%
      filter(Group != "Other") %>%
      group_by(Group) %>%
      summarise(Total = sum(Mean_Pct), .groups = "drop") %>%
      arrange(desc(Total))
    print(group_summary)
    
    # 11. Output Top variables
    cat("\n--- Top 10 most important variables ---\n")
    print(head(avg_imp, 10))
  }
  cat("\n========================================\n")
  cat(sprintf("  ✅ Complete! Analyzed %d response variables\n", length(all_results)))
  cat(sprintf("  📁 Output directory: %s\n", out_dir))
  cat("========================================\n")
  
  return(all_results)
}

# ============================================
# 7. Run all groups
# ============================================

run_all_grouped_importance_by_zn <- function() {
  
  cat("\n========================================\n")
  cat("  Running grouped variable importance analysis (three groups)\n")
  cat("  All / Low Zn / High Zn | Scheme A standard meta random forest\n")
  cat("========================================\n")
  
  # Define groups
  zn_configs <- list(
    list(level = NULL, label = "All", file = "all_processed_data.rds"),
    list(level = "low", label = "Low Zn Level", file = "all_processed_data_low.rds"),
    list(level = "high", label = "High Zn Level", file = "all_processed_data_high.rds")
  )
  
  all_group_results <- list()
  
  for(config in zn_configs) {
    level <- config$level
    label <- config$label
    file_name <- config$file
    
    cat(sprintf("\n========== Processing %s ==========\n", label))
    
    file_path <- file.path(PROCESSED_DIR, file_name)
    
    if(!file.exists(file_path)) {
      cat(sprintf("⚠️ File does not exist: %s, skipping\n", file_path))
      next
    }
    
    processed_data <- readRDS(file_path)
    cat(sprintf("✓ Loaded data: %s (%d response variables)\n", file_name, length(processed_data)))
    
    # Run analysis
    results <- run_all_grouped_importance(processed_data, level, label)
    all_group_results[[label]] <- results
  }
  
  cat("\n========================================\n")
  cat("  ✅ All group analyses complete!\n")
  cat("========================================\n")
  cat("\nOutput directories:\n")
  cat("  results/figures/grouped_importance/all/     - All group\n")
  cat("  results/figures/grouped_importance/low/     - Low Zn group\n")
  cat("  results/figures/grouped_importance/high/    - High Zn group\n")
  
  return(all_group_results)
}

# ============================================
# Run entry
# ============================================
if(USE_ZN_LEVEL_GROUP) {
  # Strongly recommended to clear environment before running to avoid old data interference
  # rm(list=ls());gc()
  grouped_results <- run_all_grouped_importance_by_zn()
} else {
  cat("USE_ZN_LEVEL_GROUP not enabled, running All group only\n")
  if(file.exists(file.path(PROCESSED_DIR, "all_processed_data.rds"))) {
    processed_data <- readRDS(file.path(PROCESSED_DIR, "all_processed_data.rds"))
    grouped_results <- run_all_grouped_importance(processed_data, NULL, "All")
  }
}