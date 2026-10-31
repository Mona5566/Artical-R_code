# ============================================
# Winsorization script: 00_winsorize_data.R (supports Zn_Level grouping)
# Function: winsorize computed lnRR effect sizes overall
#           save full data and grouped data after winsorization
# ============================================

source("00_config.R")

# ========== Corrected winsorization function (two-sided) ==========
winsorize_lnRR <- function(x, prob = 0.95) {
  # prob = 0.95 means retaining the middle 95% of data
  # i.e., winsorize to 2.5% and 97.5% quantiles
  lower_prob <- (1 - prob) / 2  # 0.025
  upper_prob <- 1 - lower_prob   # 0.975
  
  lower <- quantile(x, lower_prob, na.rm = TRUE)
  upper <- quantile(x, upper_prob, na.rm = TRUE)
  
  # Record number of adjustments
  n_lower <- sum(x < lower, na.rm = TRUE)
  n_upper <- sum(x > upper, na.rm = TRUE)
  
  # Perform winsorization
  x_wins <- x
  x_wins[!is.na(x_wins) & x_wins < lower] <- lower
  x_wins[!is.na(x_wins) & x_wins > upper] <- upper
  
  # Print information
  cat(sprintf("    Winsorization thresholds: [%.4f, %.4f] (%.1f%%-%.1f%%)\n", 
              lower, upper, lower_prob*100, upper_prob*100))
  cat(sprintf("    Adjusted: lower %d, upper %d, total %d (%.1f%%)\n",
              n_lower, n_upper, n_lower + n_upper, 
              100 * (n_lower + n_upper) / sum(!is.na(x))))
  
  return(list(
    lnRR_wins = x_wins,
    n_lower = n_lower,
    n_upper = n_upper,
    lower_thr = lower,
    upper_thr = upper
  ))
}

# ========== Execute winsorization ==========
cat("\n")
cat("========================================\n")
cat("     lnRR effect size winsorization\n")
cat("     Group by Zn_Level after winsorization\n")
cat("========================================\n")
cat(sprintf("Winsorization setting: retain middle %.0f%% of data\n", WINSORIZE_PERCENTILE * 100))
cat(sprintf("i.e., winsorize to %.1f%% and %.1f%% quantiles\n", 
            (1 - WINSORIZE_PERCENTILE)/2 * 100,
            (1 - (1 - WINSORIZE_PERCENTILE)/2) * 100))
cat(sprintf("Zn_Level grouping: %s\n", ifelse(USE_ZN_LEVEL_GROUP, "Yes (grouped after winsorization)", "No")))
cat("========================================\n\n")

# Create directory for winsorized data
WINSORIZED_DIR <- "data/winsorized/"
dir.create(WINSORIZED_DIR, recursive = TRUE, showWarnings = FALSE)

# Log for winsorization processing
winsorize_log <- data.frame(
  response = character(),
  n_total = integer(),
  n_lower = integer(),
  n_upper = integer(),
  lower_threshold = numeric(),
  upper_threshold = numeric(),
  original_min = numeric(),
  original_max = numeric(),
  winsorized_min = numeric(),
  winsorized_max = numeric(),
  stringsAsFactors = FALSE
)

# Step 1: compute effect sizes for all response variables (using raw data)
cat("Step 1: Computing lnRR effect sizes...\n\n")

# Temporarily set data source to raw data
DATA_DIR_ORIG <- DATA_DIR
DATA_DIR <- "data/raw/"

# Run effect size computation
source("01_effect_size.R", local = TRUE)

# Restore DATA_DIR
DATA_DIR <- DATA_DIR_ORIG

# Check whether processed_data exists
if(!exists("processed_data") || length(processed_data) == 0) {
  stop("Error: processed_data not generated, please check 01_effect_size.R")
}

# Step 2: winsorize lnRR for each response variable overall
cat("\n\nStep 2: Winsorizing lnRR overall...\n")

all_data_winsorized <- list()

for(i in 1:nrow(RESPONSE_FILES)) {
  resp <- RESPONSE_FILES$response_name[i]
  dat <- processed_data[[resp]]
  
  if(is.null(dat)) {
    cat(sprintf("\nSkipping: %s (data not found)\n", resp))
    next
  }
  
  cat(sprintf("\nProcessing: %s\n", resp))
  
  # Check whether Zn_Level column exists
  if(!"Zn_Level" %in% colnames(dat)) {
    cat("  Warning: data lacks Zn_Level column\n")
    # If missing, add default value
    dat$Zn_Level <- "all"
  }
  
  # Print Zn_Level distribution (before winsorization)
  zn_counts <- table(dat$Zn_Level)
  cat(sprintf("  Zn_Level distribution (before winsorization): %s\n", 
              paste(names(zn_counts), "=", zn_counts, collapse = ", ")))
  
  cat(sprintf("  Original lnRR range: [%.4f, %.4f]\n", 
              min(dat$lnRR, na.rm = TRUE), 
              max(dat$lnRR, na.rm = TRUE)))
  
  # Save original values
  original_min <- min(dat$lnRR, na.rm = TRUE)
  original_max <- max(dat$lnRR, na.rm = TRUE)
  n_total <- sum(!is.na(dat$lnRR))
  
  # Save original lnRR
  dat$lnRR_original <- dat$lnRR
  
  # ========== Overall winsorization (not by group) ==========
  cat("  Performing overall winsorization...\n")
  winsorized <- winsorize_lnRR(dat$lnRR, prob = WINSORIZE_PERCENTILE)
  dat$lnRR <- winsorized$lnRR_wins
  
  # Record log
  winsorize_log <- rbind(winsorize_log, data.frame(
    response = resp,
    n_total = n_total,
    n_lower = winsorized$n_lower,
    n_upper = winsorized$n_upper,
    lower_threshold = winsorized$lower_thr,
    upper_threshold = winsorized$upper_thr,
    original_min = original_min,
    original_max = original_max,
    winsorized_min = min(dat$lnRR, na.rm = TRUE),
    winsorized_max = max(dat$lnRR, na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
  
  # Recompute related effect sizes
  dat <- dat %>%
    mutate(
      pct_change = (exp(lnRR) - 1) * 100,
      lnRR_lower = lnRR - 1.96 * sqrt(variance),
      lnRR_upper = lnRR + 1.96 * sqrt(variance),
      pct_lower = (exp(lnRR_lower) - 1) * 100,
      pct_upper = (exp(lnRR_upper) - 1) * 100
    )
  
  cat(sprintf("  Winsorized lnRR range: [%.4f, %.4f]\n", 
              min(dat$lnRR, na.rm = TRUE), 
              max(dat$lnRR, na.rm = TRUE)))
  
  # Print Zn_Level distribution after winsorization
  zn_counts_after <- table(dat$Zn_Level)
  cat(sprintf("  Zn_Level distribution (after winsorization): %s\n", 
              paste(names(zn_counts_after), "=", zn_counts_after, collapse = ", ")))
  
  # Save winsorized full data to CSV
  output_path <- file.path(WINSORIZED_DIR, RESPONSE_FILES$file[i])
  write_csv(dat, output_path)
  cat(sprintf("  Saved full data: %s\n", output_path))
  
  all_data_winsorized[[resp]] <- dat
}

# Save winsorization log
log_path <- file.path(WINSORIZED_DIR, "winsorize_lnRR_log.csv")
write_csv(winsorize_log, log_path)
cat(sprintf("\nWinsorization log saved: %s\n", log_path))

# Print summary
cat("\n========== Winsorization Summary ==========\n")
print(winsorize_log[, c("response", "n_total", "n_lower", "n_upper", 
                        "lower_threshold", "upper_threshold")])

# ========== Step 3: save full data and grouped data ==========
cat("\n\nStep 3: Saving full data and grouped data...\n")

# 1. Save full data (after winsorization)
saveRDS(all_data_winsorized, file.path(PROCESSED_DIR, "all_processed_data.rds"))
cat(sprintf("✅ Saved full data: all_processed_data.rds (%d response variables)\n", 
            length(all_data_winsorized)))

# 2. Save grouped data (grouped by Zn_Level after winsorization)
if(USE_ZN_LEVEL_GROUP) {
  
  all_data_low <- list()
  all_data_high <- list()
  
  for(resp in names(all_data_winsorized)) {
    dat <- all_data_winsorized[[resp]]
    
    # Check whether Zn_Level column exists
    if("Zn_Level" %in% colnames(dat)) {
      # Low Zn group
      dat_low <- dat %>% filter(Zn_Level == "low")
      if(nrow(dat_low) > 0) {
        all_data_low[[resp]] <- dat_low
      }
      
      # High Zn group
      dat_high <- dat %>% filter(Zn_Level == "high")
      if(nrow(dat_high) > 0) {
        all_data_high[[resp]] <- dat_high
      }
    }
  }
  
  # Save Low Zn group data
  if(length(all_data_low) > 0) {
    saveRDS(all_data_low, file.path(PROCESSED_DIR, "all_processed_data_low.rds"))
    cat(sprintf("✅ Saved Low Zn group data: all_processed_data_low.rds (%d response variables)\n", 
                length(all_data_low)))
    
    # Print Low Zn sample size for each response variable
    low_n <- sapply(all_data_low, nrow)
    cat("  Low Zn sample sizes:\n")
    for(resp in names(low_n)) {
      cat(sprintf("    %s: %d\n", resp, low_n[resp]))
    }
  } else {
    cat("⚠️ No Low Zn group data found\n")
  }
  
  # Save High Zn group data
  if(length(all_data_high) > 0) {
    saveRDS(all_data_high, file.path(PROCESSED_DIR, "all_processed_data_high.rds"))
    cat(sprintf("✅ Saved High Zn group data: all_processed_data_high.rds (%d response variables)\n", 
                length(all_data_high)))
    
    # Print High Zn sample size for each response variable
    high_n <- sapply(all_data_high, nrow)
    cat("  High Zn sample sizes:\n")
    for(resp in names(high_n)) {
      cat(sprintf("    %s: %d\n", resp, high_n[resp]))
    }
  } else {
    cat("⚠️ No High Zn group data found\n")
  }
  
} else {
  cat("⚠️ USE_ZN_LEVEL_GROUP = FALSE, skipping grouped data saving\n")
}

cat("\n========== lnRR Winsorization Complete ==========\n")
cat("\nGenerated output files:\n")
cat("  data/winsorized/*.csv           - Winsorized CSV files\n")
cat("  data/processed/all_processed_data.rds      - Full data\n")
if(USE_ZN_LEVEL_GROUP) {
  cat("  data/processed/all_processed_data_low.rds  - Low Zn group data\n")
  cat("  data/processed/all_processed_data_high.rds - High Zn group data\n")
}
cat("  data/winsorized/winsorize_lnRR_log.csv - Winsorization log\n")

cat("\nYou can now run modules 2-10 for analysis\n")