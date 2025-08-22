################################################################################
#### RESULTS FOR TABLES ########################################################
################################################################################

library(dplyr); library(tidyr); library(stringr); library(tibble)

# thresholds used
PSR_BAD_TH <- 1.05
ESS_LOW_TH <- 400
EPS_TRUTH  <- 1e-12 # to avoid division by 0 (rel bias)
DIF_SIZE   <- 0.30
ALPHA_P    <- 0.05

.get_param_names <- function(df) {
  rn <- rownames(df)
  if (!is.null(rn) && any(nchar(rn) > 0)) return(rn)
  nm <- tolower(names(df))
  pick <- match(TRUE, nm %in% c("param","parameter","name","label"), nomatch = 0L)
  if (pick > 0) return(as.character(df[[pick]]))
  rep(NA_character_, nrow(df))
}

.extract_item <- function(x) {
  m <- stringr::str_match(tolower(x), "y([2-5])")[,2]
  ifelse(is.na(m), NA_character_, paste0("y", m))
}

# one interaction row per item using (explicit b_* or non-standardized "y ~ X*lat")
.match_interaction_rows <- function(param_names, family) {
  if (is.null(param_names)) return(integer(0))
  p <- tolower(param_names)
  items <- paste0("y", 2:5)
  idx_keep <- integer(0)
  
  for (it in items) {
    hits <- if (identical(family, "latX")) {
      h <- grep(paste0("\\bb_int_", it, "\\b"), p, perl = TRUE)                 # explicit
      if (!length(h)) h <- grep(paste0("^", it, "\\s*~\\s*latx\\*lat\\b"), p, perl = TRUE)  # label
      h
    } else {
      h <- grep(paste0("\\bb_gint_", it, "\\b"), p, perl = TRUE)                # explicit
      if (!length(h)) h <- grep(paste0("^", it, "\\s*~\\s*group\\*lat\\b"), p, perl = TRUE) # label
      h
    }
    if (length(hits)) idx_keep <- c(idx_keep, hits[1])
  }
  unique(idx_keep)
}

.true_val <- function(item, dif_items, family, dif_size = DIF_SIZE) {
  if (identical(family, "latX") && item %in% dif_items) dif_size else 0.0
}

# core results extraction 
make_est_rows <- function(node, family, prior_id, meta, dif_items) {
  est_raw <- node$estimates_raw
  if (is.null(est_raw)) return(NULL)
  
  df <- as.data.frame(est_raw, check.names = FALSE)  # keep "2.5%"/"97.5%"
  pnames <- .get_param_names(df)
  idx <- .match_interaction_rows(pnames, family)
  if (!length(idx)) return(NULL)
  
  df_sub <- df[idx, , drop = FALSE]
  items  <- .extract_item(pnames[idx])
  
  tibble(
    N = meta$N,
    dif_case = if (identical(sort(dif_items), c("y2","y3"))) "2DIF" else "AllDIF",
    model_family = family,
    prior_id = prior_id,
    rep = meta$rep,
    item = items,
    
    # exact column names as from rBlimp (mapping didnt work well)
    Mean     = df_sub$Mean,
    Median   = df_sub$Median,
    StdDev   = df_sub$StdDev,
    `2.5%`   = df_sub$`2.5%`,
    `97.5%`  = df_sub$`97.5%`,
    PSR      = df_sub$PSR,
    N_Eff    = df_sub$N_Eff,
    MCMC_SE  = df_sub$MCMC_SE,
    ChiSq    = if ("ChiSq"  %in% names(df_sub))  df_sub$ChiSq  else NA_real_,
    PValue   = suppressWarnings(as.numeric(df_sub$PValue)),
    
    truth       = vapply(items, \(it) .true_val(it, dif_items, family), numeric(1)),
    is_dif_scen = if (identical(family, "latX")) items %in% dif_items else FALSE
  )
}

collect_est_rows <- function(results_sim) {
  out <- list()
  for (scen_key in names(results_sim)) {
    reps <- results_sim[[scen_key]]
    for (r in seq_along(reps)) {
      meta <- reps[[r]]$meta; meta$rep <- reps[[r]]$rep
      dif_items <- meta$dif_items
      
      if (!is.null(reps[[r]]$latX$full)) {
        for (prior_id in names(reps[[r]]$latX$full)) {
          out[[length(out)+1]] <- make_est_rows(
            reps[[r]]$latX$full[[prior_id]], "latX", prior_id, meta, dif_items
          )
        }
      }
      if (!is.null(reps[[r]]$group$full)) {
        for (prior_id in names(reps[[r]]$group$full)) {
          out[[length(out)+1]] <- make_est_rows(
            reps[[r]]$group$full[[prior_id]], "group", prior_id, meta, dif_items
          )
        }
      }
    }
  }
  dplyr::bind_rows(out)
}

est_wide <- collect_est_rows(results_sim)

# long for mean and median (intervals/p-values are the same per param, but we keep both estimators for completeness)
est_long <- est_wide %>%
  tidyr::pivot_longer(c("Mean","Median"), names_to = "estimator", values_to = "estimate")

perf_min <- est_long %>%
  group_by(N, dif_case, model_family, prior_id, item, estimator) %>%
  summarise(
    reps    = n_distinct(rep),
    is_dif  = dplyr::first(is_dif_scen),
    
    # accuracy
    bias     = mean(estimate - truth, na.rm = TRUE),
    rel_bias = if (abs(first(truth)) > EPS_TRUTH)
      mean((estimate - truth)/abs(truth), na.rm = TRUE) else NA_real_,
    
    # coverage
    coverage95 = mean(`2.5%` <= truth & `97.5%` >= truth, na.rm = TRUE),
    
    # SD/SE calibration
    emp_sd       = sd(estimate, na.rm = TRUE),
    StdDev_mean  = mean(StdDev, na.rm = TRUE),
    sd_over_se   = emp_sd / pmax(StdDev_mean, 1e-12),
    
    # --- detection (both ways) ---
    # CI-based: does the 95% CI exclude 0?
    power_ci       = mean(ifelse(first(is_dif_scen), (`2.5%` > 0 | `97.5%` < 0), NA), na.rm = TRUE),
    type1_error_ci = mean(ifelse(!first(is_dif_scen), (`2.5%` > 0 | `97.5%` < 0), NA), na.rm = TRUE),
    
    # p-value–based: is PValue < alpha?
    power_p        = mean(ifelse(first(is_dif_scen), PValue < ALPHA_P, NA), na.rm = TRUE),
    type1_error_p  = mean(ifelse(!first(is_dif_scen), PValue < ALPHA_P, NA), na.rm = TRUE),
    
    # diagnostics
    psr_mean    = mean(PSR,     na.rm = TRUE),
    psr_median  = median(PSR,   na.rm = TRUE),
    ess_mean    = mean(N_Eff,   na.rm = TRUE),
    ess_median  = median(N_Eff, na.rm = TRUE),
    mcse_mean   = mean(MCMC_SE, na.rm = TRUE),
    
    psr_bad_prop_105 = mean(PSR > PSR_BAD_TH,   na.rm = TRUE),
    ess_low_prop_400 = mean(N_Eff < ESS_LOW_TH, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(N, dif_case, model_family, prior_id, item, estimator)
