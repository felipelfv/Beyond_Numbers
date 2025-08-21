# --- Minimal PSR checker -----------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tibble) })
stopifnot(exists("results_sim"))

check_psr_convergence <- function(results_sim, threshold = 1.05) {
  stopifnot(!is.null(results_sim))
  
  # normalize PSR column name; prefer your .normalize_cols() if present
  norm_cols <- if (exists(".normalize_cols")) .normalize_cols else function(df) {
    if (is.null(df) || !is.data.frame(df)) return(df)
    nm <- tolower(names(df))
    nm[nm %in% c("rhat","psrf")] <- "psr"
    names(df) <- nm
    df
  }
  
  add_fit_row <- function(node, meta, fam, model, prior_id = NA_character_, threshold = 1.05) {
    df <- try(as.data.frame(node$estimates_raw, check.names = FALSE), silent = TRUE)
    if (inherits(df, "try-error") || is.null(df) || !nrow(df)) return(NULL)
    df <- norm_cols(df)
    if (!("psr" %in% names(df))) return(NULL)
    
    psr <- suppressWarnings(as.numeric(df$psr))
    n_total <- sum(!is.na(psr))
    if (n_total == 0) return(NULL)
    
    tibble(
      N              = meta$N,
      dif_case       = if (identical(sort(meta$dif_items), c("y2","y3"))) "2DIF" else "AllDIF",
      model_family   = fam,
      model          = model,                  # "full" or "null"
      prior_id       = prior_id,
      rep            = if (!is.null(meta$rep)) meta$rep else NA_integer_,
      n_params_psr   = n_total,
      n_bad_psr      = sum(psr > threshold, na.rm = TRUE),
      max_psr        = max(psr, na.rm = TRUE),
      converged      = max_psr <= threshold
    )
  }
  
  rows <- list()
  for (scen in names(results_sim)) {
    reps <- results_sim[[scen]]
    for (r in seq_along(reps)) {
      meta <- reps[[r]]$meta; meta$rep <- reps[[r]]$rep
      for (fam in c("latX","group")) {
        # FULL models (all priors)
        full_priors <- names(reps[[r]][[fam]]$full)
        for (pid in full_priors) {
          rows[[length(rows)+1]] <- add_fit_row(reps[[r]][[fam]]$full[[pid]], meta, fam, "full", pid, threshold)
        }
        # NULL model
        rows[[length(rows)+1]] <- add_fit_row(reps[[r]][[fam]]$null, meta, fam, "null", NA_character_, threshold)
      }
    }
  }
  
  by_fit <- dplyr::bind_rows(rows)
  
  rates <- by_fit %>%
    group_by(N, dif_case, model_family, model, prior_id) %>%
    summarise(
      reps           = n(),
      conv_rate      = mean(converged, na.rm = TRUE),     # proportion with max PSR <= threshold
      any_psr_gt_thr = mean(n_bad_psr > 0, na.rm = TRUE), # ~= 1 - conv_rate when no NAs
      .groups = "drop"
    )
  
  list(by_fit = by_fit, rates = rates)
}

# --- Run + quick views -------------------------------------------------------
psr_chk <- check_psr_convergence(results_sim, threshold = 1.05)

# per-fit details (one row per rep × family × model × prior)
dplyr::glimpse(psr_chk$by_fit)

# convergence rates (paper-style)
psr_chk$rates %>% arrange(desc(conv_rate))

# overall yes/no
all_ok <- all(psr_chk$by_fit$converged)
all_ok

# quick counts
dplyr::count(psr_chk$by_fit, converged)

# failures only (if any)
failures <- psr_chk$by_fit %>%
  dplyr::filter(!converged) %>%
  dplyr::arrange(dplyr::desc(max_psr))
failures

if (all(psr_chk$by_fit$converged)) {
  message("All fits converged (max PSR \u2264 1.05).")
} else {
  print(failures %>% dplyr::select(N, dif_case, model_family, model, prior_id, rep, n_bad_psr, max_psr))
}
