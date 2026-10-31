# ============================================
# Fix Meta regression table generation issue
# ============================================

library(tidyverse)
library(metafor)
library(patchwork)
library(viridis)

# Set paths
if(file.exists("00_config.R")) {
  source("00_config.R")
} else {
  PROCESSED_DIR <- "data/processed"
  FIGURES_DIR <- "figures"
  TABLES_DIR <- "tables"
  dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
}

# Load data
data_file <- file.path(PROCESSED_DIR, "all_processed_data.rds")
if(file.exists(data_file)) {
  processed_data <- readRDS(data_file)
  cat("✅ Data loaded successfully!\n")
} else {
  stop("❌ Data file does not exist: ", data_file)
}

# ============================================
# fit_model function
# ============================================
fit_model <- function(formula_mods, data, use_random = TRUE) {
  if(use_random && length(unique(data$study_id)) >= 3) {
    tryCatch({
      return(rma.mv(
        yi = lnRR, V = variance,
        random = ~ 1 | study_id,
        mods = formula_mods,
        data = data,
        method = "REML",
        control = list(optimizer = "optim", optmethod = "Nelder-Mead",
                       maxit = 2000)
      ))
    }, error = function(e) {
      tryCatch({
        return(rma(
          yi = lnRR, vi = variance,
          mods = formula_mods,
          data = data,
          method = "REML"
        ))
      }, error = function(e2) {
        return(NULL)
      })
    })
  } else {
    tryCatch({
      return(rma(
        yi = lnRR, vi = variance,
        mods = formula_mods,
        data = data,
        method = "REML"
      ))
    }, error = function(e) {
      return(NULL)
    })
  }
}

# ============================================
# compare_metareg_models function (simplified version)
# ============================================
compare_metareg_models <- function(dat, moderator, response_name) {
  
  # Data preparation
  if(moderator == "Total_Zn") {
    if(any(dat[[moderator]] <= 0, na.rm = TRUE)) {
      dat$mod_transformed <- log10(dat[[moderator]] + 0.01)
    } else {
      dat$mod_transformed <- log10(dat[[moderator]])
    }
    x_label <- expression("Total Soil Zn (mg kg"^{-1} * ", log"[10]~"scale)")
  } else {
    dat$mod_transformed <- dat[[moderator]]
    x_label <- moderator
  }
  
  # Clean data
  dat_clean <- dat %>%
    filter(
      !is.na(mod_transformed),
      !is.na(lnRR),
      is.finite(lnRR),
      variance > 0,
      is.finite(variance),
      !is.na(study_id)
    ) %>%
    group_by(study_id) %>%
    filter(n() >= 2) %>%
    ungroup()
  
  if(nrow(dat_clean) < 10) {
    return(NULL)
  }
  
  if(length(unique(dat_clean$mod_transformed)) < 3) {
    return(NULL)
  }
  
  n_studies <- length(unique(dat_clean$study_id))
  
  # Linear model
  model_linear <- fit_model(~ mod_transformed, dat_clean, use_random = n_studies >= 3)
  if(is.null(model_linear)) {
    return(NULL)
  }
  
  aic_linear <- tryCatch(AIC(model_linear), error = function(e) NA)
  p_linear <- tryCatch(model_linear$pval[2], error = function(e) NA)
  linear_valid <- !is.na(p_linear) && p_linear < 0.05
  
  # Quadratic model
  model_quad <- NULL
  aic_quad <- NA
  p_quad <- NA
  quad_valid <- FALSE
  vertex <- NA
  vertex_type <- NA
  
  if(nrow(dat_clean) >= 20) {
    x_mean <- mean(dat_clean$mod_transformed, na.rm = TRUE)
    dat_clean$mod_centered <- dat_clean$mod_transformed - x_mean
    dat_clean$mod_sq <- dat_clean$mod_centered^2
    
    cor_check <- tryCatch(
      cor(dat_clean$mod_centered, dat_clean$mod_sq, use = "complete.obs"),
      error = function(e) NA
    )
    
    if(!is.na(cor_check) && abs(cor_check) < 0.999) {
      model_quad <- fit_model(~ mod_centered + mod_sq, dat_clean, use_random = n_studies >= 3)
      
      if(!is.null(model_quad)) {
        aic_quad <- tryCatch(AIC(model_quad), error = function(e) NA)
        p_quad <- tryCatch(model_quad$pval[3], error = function(e) NA)
        quad_valid <- !is.na(p_quad) && p_quad < 0.05
        
        a <- tryCatch(model_quad$b[3], error = function(e) NA)
        b <- tryCatch(model_quad$b[2], error = function(e) NA)
        if(!is.na(a) && !is.na(b) && a != 0) {
          vertex <- -b / (2 * a) + x_mean
          vertex_type <- ifelse(a > 0, "minimum", "maximum")
        }
      }
    }
  }
  
  # Log model
  model_log <- NULL
  aic_log <- NA
  p_log <- NA
  log_valid <- FALSE
  
  if(moderator == "Total_Zn" && all(dat_clean$mod_transformed > 0, na.rm = TRUE)) {
    dat_clean$mod_log <- log(dat_clean$mod_transformed + 1)
    model_log <- fit_model(~ mod_log, dat_clean, use_random = n_studies >= 3)
    
    if(!is.null(model_log)) {
      aic_log <- tryCatch(AIC(model_log), error = function(e) NA)
      p_log <- tryCatch(model_log$pval[2], error = function(e) NA)
      log_valid <- !is.na(p_log) && p_log < 0.05
    }
  }
  
  # Segmented model
  model_seg <- NULL
  aic_seg <- NA
  p_seg <- NA
  seg_valid <- FALSE
  breakpoint <- NA
  
  if(nrow(dat_clean) >= 30) {
    bp_candidates <- tryCatch(
      quantile(dat_clean$mod_transformed, probs = seq(0.25, 0.75, 0.1)),
      error = function(e) NULL
    )
    
    if(!is.null(bp_candidates) && length(bp_candidates) > 0) {
      best_model <- NULL
      best_aic <- Inf
      best_bp <- NA
      best_p <- 1
      
      for(bp in bp_candidates) {
        dat_clean$seg_high <- as.numeric(dat_clean$mod_transformed > bp)
        dat_clean$mod_seg <- dat_clean$mod_transformed * dat_clean$seg_high
        
        cor_check <- tryCatch(
          cor(dat_clean$mod_transformed, dat_clean$mod_seg, use = "complete.obs"),
          error = function(e) NA
        )
        
        if(!is.na(cor_check) && abs(cor_check) > 0.99) next
        
        temp <- fit_model(~ mod_transformed + mod_seg, dat_clean, use_random = n_studies >= 3)
        
        if(!is.null(temp) && length(temp$pval) >= 3) {
          p_int <- tryCatch(temp$pval[3], error = function(e) NA)
          if(!is.na(p_int) && p_int < 0.05) {
            aic_temp <- tryCatch(AIC(temp), error = function(e) Inf)
            if(!is.infinite(aic_temp) && aic_temp < best_aic) {
              best_aic <- aic_temp
              best_model <- temp
              best_bp <- bp
              best_p <- p_int
            }
          }
        }
      }
      
      if(!is.null(best_model)) {
        model_seg <- best_model
        aic_seg <- best_aic
        p_seg <- best_p
        breakpoint <- best_bp
        seg_valid <- TRUE
      }
    }
  }
  
  # Model selection
  models_list <- list()
  
  if(linear_valid && !is.null(model_linear)) {
    models_list[["Linear"]] <- list(
      model = model_linear, aic = aic_linear, p = p_linear, type = "Linear"
    )
  }
  if(quad_valid && !is.null(model_quad)) {
    models_list[["Quadratic"]] <- list(
      model = model_quad, aic = aic_quad, p = p_quad, type = "Quadratic",
      vertex = vertex, vertex_type = vertex_type
    )
  }
  if(log_valid && !is.null(model_log)) {
    models_list[["Log"]] <- list(
      model = model_log, aic = aic_log, p = p_log, type = "Log"
    )
  }
  if(seg_valid && !is.null(model_seg)) {
    models_list[["Segmented"]] <- list(
      model = model_seg, aic = aic_seg, p = p_seg, type = "Segmented",
      breakpoint = breakpoint
    )
  }
  
  # If no significant model, use linear model (even if not significant)
  if(length(models_list) == 0) {
    if(!is.null(model_linear)) {
      selected_model <- model_linear
      model_type <- "Linear (non-significant)"
      best_name <- "Linear_NS"
    } else {
      return(NULL)
    }
  } else {
    best_name <- names(models_list)[which.min(sapply(models_list, function(x) x$aic))]
    best <- models_list[[best_name]]
    selected_model <- best$model
    model_type <- best$type
  }
  
  # Return results
  return(list(
    selected_model = selected_model,
    model_type = model_type,
    best_name = if(exists("best_name")) best_name else "Linear_NS",
    all_models = models_list,
    dat_clean = dat_clean,
    moderator = moderator,
    response_name = response_name,
    x_label = x_label,
    n = nrow(dat_clean),
    n_studies = n_studies,
    linear = list(p = if(exists("p_linear")) p_linear else NA, 
                  aic = if(exists("aic_linear")) aic_linear else NA),
    quadratic = list(p = if(exists("p_quad")) p_quad else NA, 
                     aic = if(exists("aic_quad")) aic_quad else NA, 
                     valid = if(exists("quad_valid")) quad_valid else FALSE,
                     vertex = if(exists("vertex")) vertex else NA, 
                     vertex_type = if(exists("vertex_type")) vertex_type else NA),
    segmented = list(p = if(exists("p_seg")) p_seg else NA, 
                     aic = if(exists("aic_seg")) aic_seg else NA, 
                     valid = if(exists("seg_valid")) seg_valid else FALSE,
                     breakpoint = if(exists("breakpoint")) breakpoint else NA),
    log = list(p = if(exists("p_log")) p_log else NA, 
               aic = if(exists("aic_log")) aic_log else NA, 
               valid = if(exists("log_valid")) log_valid else FALSE)
  ))
}

# ============================================
# Run analysis and generate tables
# ============================================

cat("\n", rep("=", 70), "\n", sep="")
cat("🧠 Smart Meta-regression - Table fixed version\n")
cat(rep("=", 70), "\n", sep="")

responses <- c(
  "Zn_shoot", "Zn_root", "Zn_plant",
  "Biomass_shoot", "Biomass_root", "Biomass_fruit", 
  "Biomass_leave", "Plant_height",
  "BCF_shoot", "BCF_root", "BCF_plant",
  "P_shoot", "P_root"
)

moderators <- c("Total_Zn", "soil_pH")

all_results <- list()
success_count <- list(Total_Zn = 0, soil_pH = 0)
total_count <- list(Total_Zn = 0, soil_pH = 0)

for(resp in responses) {
  if(!resp %in% names(processed_data)) next
  
  dat <- processed_data[[resp]]
  if(is.null(dat)) next
  
  for(mod in moderators) {
    if(!mod %in% colnames(dat)) next
    
    total_count[[mod]] <- total_count[[mod]] + 1
    
    cat(sprintf("\n📊 %s ~ %s", resp, mod))
    
    model_fits <- compare_metareg_models(dat, mod, resp)
    
    if(is.null(model_fits)) {
      cat(" ❌ Failed\n")
      next
    }
    
    success_count[[mod]] <- success_count[[mod]] + 1
    all_results[[paste(resp, mod, sep = "_")]] <- model_fits
    cat(sprintf(" ✅ (n=%d, studies=%d, model=%s)", 
                model_fits$n, model_fits$n_studies, model_fits$model_type))
  }
}

cat("\n\n", rep("=", 70), "\n", sep="")
cat("📊 Generating summary table\n")
cat(rep("=", 70), "\n", sep="")

# ============================================
# Generate summary table (fixed version)
# ============================================

if(length(all_results) > 0) {
  
  cat(sprintf("Total %d successful results\n", length(all_results)))
  
  # Build table using loop
  summary_rows <- list()
  row_index <- 1
  
  for(name in names(all_results)) {
    res <- all_results[[name]]
    
    # Parse name
    parts <- strsplit(name, "_")[[1]]
    resp <- parts[1]
    mod <- paste(parts[-1], collapse = "_")
    
    cat(sprintf("\nProcessing: %s", name))
    
    # Safe extraction function
    safe_extract <- function(x, default = NA) {
      tryCatch(x, error = function(e) default)
    }
    
    # Extract p-values and AIC for all models
    linear_p <- safe_extract(res$linear$p)
    linear_aic <- safe_extract(res$linear$aic)
    quad_p <- safe_extract(res$quadratic$p)
    quad_aic <- safe_extract(res$quadratic$aic)
    seg_p <- safe_extract(res$segmented$p)
    seg_aic <- safe_extract(res$segmented$aic)
    log_p <- safe_extract(res$log$p)
    log_aic <- safe_extract(res$log$aic)
    
    # Extract vertex and breakpoint
    vertex <- safe_extract(res$quadratic$vertex)
    vertex_type <- safe_extract(res$quadratic$vertex_type)
    breakpoint <- safe_extract(res$segmented$breakpoint)
    
    # Extract validity
    quad_valid <- safe_extract(res$quadratic$valid, FALSE)
    seg_valid <- safe_extract(res$segmented$valid, FALSE)
    log_valid <- safe_extract(res$log$valid, FALSE)
    
    # Compute minimum p-value (from all models)
    all_p <- c(linear_p, quad_p, seg_p, log_p)
    all_p <- all_p[!is.na(all_p) & is.finite(all_p)]
    min_p <- if(length(all_p) > 0) min(all_p) else NA
    
    # Extract best model AIC
    best_aic <- tryCatch(AIC(res$selected_model), error = function(e) NA)
    
    # Extract slope
    slope <- tryCatch({
      coef(res$selected_model)[2]
    }, error = function(e) NA)
    
    # Build row
    row <- tryCatch({
      tibble(
        Response = resp,
        Moderator = mod,
        N = res$n,
        N_Studies = res$n_studies,
        Best_Model = res$model_type,
        Best_AIC = round(best_aic, 1),
        Best_p = round(min_p, 4),
        Slope = round(slope, 4),
        Linear_p = round(linear_p, 4),
        Linear_AIC = round(linear_aic, 1),
        Quadratic_p = round(quad_p, 4),
        Quadratic_AIC = round(quad_aic, 1),
        Quadratic_Valid = quad_valid,
        Vertex = round(vertex, 2),
        Vertex_Type = ifelse(!is.na(vertex_type), vertex_type, NA),
        Segmented_p = round(seg_p, 4),
        Segmented_AIC = round(seg_aic, 1),
        Segmented_Valid = seg_valid,
        Breakpoint = round(breakpoint, 2),
        Log_p = round(log_p, 4),
        Log_AIC = round(log_aic, 1),
        Log_Valid = log_valid
      )
    }, error = function(e) {
      cat(sprintf(" ⚠️ Error: %s", e$message))
      return(NULL)
    })
    
    if(!is.null(row)) {
      summary_rows[[row_index]] <- row
      row_index <- row_index + 1
      cat(" ✅")
    }
  }
  
  # Merge all rows
  if(length(summary_rows) > 0) {
    summary_table <- bind_rows(summary_rows)
    
    cat("\n\n", rep("=", 70), "\n", sep="")
    cat("✅ Summary table generated successfully!\n")
    cat(rep("=", 70), "\n", sep="")
    
    # Save table
    tables_dir <- file.path(TABLES_DIR, "smart_metareg")
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
    write_csv(summary_table, file.path(tables_dir, "smart_metareg_results.csv"))
    cat(sprintf("\n💾 Table saved: %s\n", 
                file.path(tables_dir, "smart_metareg_results.csv")))
    
    # Print summary
    cat("\n📊 Model selection summary:\n")
    print(table(summary_table$Best_Model))
    
    cat("\n📊 Table preview:\n")
    print(summary_table %>% select(Response, Moderator, N, Best_Model, Best_AIC, Best_p))
    
    # Statistics
    cat("\n📈 Successful analysis count:\n")
    cat(sprintf("  Total_Zn: %d / %d\n", success_count$Total_Zn, total_count$Total_Zn))
    cat(sprintf("  soil_pH: %d / %d\n", success_count$soil_pH, total_count$soil_pH))
    
    # Significant model statistics
    cat("\n📊 Significant model distribution:\n")
    sig_models <- summary_table %>%
      filter(!grepl("non-significant", Best_Model, ignore.case = TRUE))
    cat(sprintf("  Significant models: %d / %d (%.1f%%)\n", 
                nrow(sig_models), nrow(summary_table),
                nrow(sig_models)/nrow(summary_table)*100))
    
    # Breakpoint statistics
    if(any(!is.na(summary_table$Breakpoint))) {
      cat("\n📍 Breakpoint statistics (segmented model):\n")
      bp_summary <- summary_table %>%
        filter(!is.na(Breakpoint)) %>%
        summarise(
          Mean = mean(Breakpoint, na.rm = TRUE),
          SD = sd(Breakpoint, na.rm = TRUE),
          Min = min(Breakpoint, na.rm = TRUE),
          Max = max(Breakpoint, na.rm = TRUE),
          n = n()
        )
      print(bp_summary)
    }
    
    # Vertex statistics
    if(any(!is.na(summary_table$Vertex))) {
      cat("\n📍 Vertex statistics (quadratic model):\n")
      v_summary <- summary_table %>%
        filter(!is.na(Vertex)) %>%
        summarise(
          Mean = mean(Vertex, na.rm = TRUE),
          SD = sd(Vertex, na.rm = TRUE),
          Min = min(Vertex, na.rm = TRUE),
          Max = max(Vertex, na.rm = TRUE),
          n = n()
        )
      print(v_summary)
    }
    
    # Save to global environment
    assign("smart_metareg_summary", summary_table, envir = .GlobalEnv)
    assign("smart_metareg_results", all_results, envir = .GlobalEnv)
    
  } else {
    cat("\n❌ No rows successfully generated\n")
  }
} else {
  cat("\n❌ No successful results\n")
}

cat("\n", rep("=", 70), "\n", sep="")
cat("✅ Analysis complete!\n")
cat(rep("=", 70), "\n", sep="")

# View generated table
if(exists("summary_table") && !is.null(summary_table)) {
  cat("\n📊 Full table:\n")
  print(summary_table)
}