# ============================================
# Module 1: Effect size calculation (supports Zn_Level grouping)
# ============================================

source("00_config.R")

compute_effect_sizes <- function(filepath, output_path = NULL) {
  
  dat <- read_csv(filepath, show_col_types = FALSE)
  
  required_cols <- c("study_id", "Treatment_n", "Control_n", 
                     "Treatment_mean", "Control_mean", 
                     "Treatment_sd", "Control_sd")
  
  missing_cols <- required_cols[!required_cols %in% colnames(dat)]
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Ensure Zn_Level exists
  if(!"Zn_Level" %in% colnames(dat)) {
    stop("Data lacks Zn_Level column, please ensure low/high group info is included")
  }
  
  # Create obs_id
  dat <- dat %>%
    arrange(study_id) %>%
    group_by(study_id) %>%
    mutate(obs_id = row_number()) %>%
    ungroup()
  
  # Compute effect sizes
  dat <- dat %>%
    mutate(
      lnRR = log(Treatment_mean / Control_mean),
      variance = (Treatment_sd^2 / (Treatment_n * Treatment_mean^2)) +
        (Control_sd^2 / (Control_n * Control_mean^2)),
      weight = 1 / variance,
      pct_change = (exp(lnRR) - 1) * 100,
      lnRR_lower = lnRR - 1.96 * sqrt(variance),
      lnRR_upper = lnRR + 1.96 * sqrt(variance),
      pct_lower = (exp(lnRR_lower) - 1) * 100,
      pct_upper = (exp(lnRR_upper) - 1) * 100
    ) %>%
    filter(is.finite(lnRR), variance > 0, variance < Inf) %>%
    mutate(
      study_id = as.factor(study_id),
      obs_id = as.factor(obs_id),
      Zn_Level = as.factor(Zn_Level)  # Ensure it is a factor
    )
  
  # Handle categorical variables
  for(mod in CATEGORICAL_MODERATORS) {
    if(mod %in% colnames(dat)) {
      dat <- dat %>% mutate(!!sym(mod) := as.factor(!!sym(mod)))
    }
  }
  
  # Handle continuous variables
  for(mod in CONTINUOUS_MODERATORS) {
    if(mod %in% colnames(dat)) {
      dat <- dat %>% mutate(!!sym(mod) := as.numeric(!!sym(mod)))
    }
  }
  
  if(is.null(output_path)) {
    output_path <- file.path(PROCESSED_DIR, basename(filepath))
  }
  write_csv(dat, output_path)
  
  # Print Zn_Level grouping info
  zn_counts <- table(dat$Zn_Level)
  cat(sprintf("✓ %s: %d observations, %d studies\n", 
              basename(filepath), nrow(dat), n_distinct(dat$study_id)))
  cat(sprintf("  Zn_Level: low=%d, high=%d\n", zn_counts["low"], zn_counts["high"]))
  
  return(dat)
}

process_all_responses <- function() {
  
  cat("\n========================================\n")
  cat("  Module 1: Effect size calculation\n")
  cat("========================================\n")
  cat(sprintf("Data source: %s\n", DATA_DIR))
  cat(sprintf("Zn_Level grouping: %s\n", ifelse(USE_ZN_LEVEL_GROUP, "Yes", "No")))
  cat("\n")
  
  all_data <- list()
  
  for(i in 1:nrow(RESPONSE_FILES)) {
    filepath <- file.path(DATA_DIR, RESPONSE_FILES$file[i])
    
    if(!file.exists(filepath)) {
      cat(sprintf("⚠ File does not exist: %s\n", filepath))
      next
    }
    
    dat <- compute_effect_sizes(filepath)
    all_data[[RESPONSE_FILES$response_name[i]]] <- dat
  }
  
  saveRDS(all_data, file.path(PROCESSED_DIR, "all_processed_data.rds"))
  
  cat("\n========================================\n")
  cat("  Data loading complete\n")
  cat("========================================\n")
  for(resp in names(all_data)) {
    df <- all_data[[resp]]
    zn_counts <- table(df$Zn_Level)
    cat(sprintf("%s: %d observations, %d studies (low=%d, high=%d)\n",
                resp, nrow(df), n_distinct(df$study_id),
                zn_counts["low"], zn_counts["high"]))
  }
  
  return(all_data)
}

# Get grouped data
get_grouped_data <- function(processed_data, group = NULL) {
  # group: NULL returns all data, "low" or "high" returns corresponding group
  
  if(is.null(group)) {
    return(processed_data)
  }
  
  grouped_data <- list()
  for(resp in names(processed_data)) {
    dat <- processed_data[[resp]]
    dat_group <- dat %>% filter(Zn_Level == group)
    if(nrow(dat_group) > 0) {
      grouped_data[[resp]] <- dat_group
    }
  }
  return(grouped_data)
}

if(interactive() || !exists("processed_data")) {
  processed_data <- process_all_responses()
  # Also save grouped data
  if(USE_ZN_LEVEL_GROUP) {
    processed_data_low <- get_grouped_data(processed_data, "low")
    processed_data_high <- get_grouped_data(processed_data, "high")
    saveRDS(processed_data_low, file.path(PROCESSED_DIR, "all_processed_data_low.rds"))
    saveRDS(processed_data_high, file.path(PROCESSED_DIR, "all_processed_data_high.rds"))
    cat("\n✅ Saved grouped data: low and high\n")
  }
}