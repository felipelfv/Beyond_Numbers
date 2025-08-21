# ==========================================================
# Minimal Monte Carlo (in-memory, simplified null storage)
# - Each replication is a list element
# - Full models keep: estimates_raw, wald_raw, modelfit_text, (optional) fit
# - Null models keep: estimates_raw, modelfit_text, (optional) fit  [no wald_raw, no "shared"]
# ==========================================================
suppressPackageStartupMessages({
  library(rblimp)
})

set.seed(12345)

# -------- CONFIG --------
priors_to_run <- c("tight","baseline","loose")
Ns    <- c(200, 500, 800)
DIF_sets <- list(`2DIF` = c("y2","y3"),
                 `AllDIF` = c("y2","y3","y4","y5"))
Y_SD <- 0.65   
n_reps <- 100

burn <- 25000
iter <- 25000
OUTPUT_SPEC <- "mean median stddev mad_sd quant95 psr n_eff mcmc_se wald pvalue"

# Store full fit objects? (large; FALSE is safer)
keep_fit_objects <- FALSE

# -------- data generator --------
gen_data <- function(N,
                     dif_items = c("y2","y3","y4","y5"),
                     p_group = 0.5, mu_shift = 0.5, sd_age = 1.0,
                     base_loading = 1.00, gamma = 0.30,
                     y_sd = 0.65,                    # <- scalar or length-5 vector
                     y_sd_range = c(0.50, 0.80)) {   # <- proper min,max default
  Group <- rbinom(N, 1, p_group)
  X     <- rnorm(N, mean = mu_shift * Group, sd = sd_age)
  Lat   <- rnorm(N, 0, 1)
  make_x <- function() 1.00*X + rnorm(N, 0, 0.55)
  x1 <- make_x(); x2 <- make_x(); x3 <- make_x(); x4 <- make_x(); x5 <- make_x()
  
  # ----- FIXED: handle y_sd correctly -----
  if (!is.null(y_sd)) {
    if (length(y_sd) == 1L) {
      sds <- rep(y_sd, 5)
    } else if (length(y_sd) == 5L) {
      sds <- y_sd
    } else {
      stop("y_sd must be NULL, a scalar, or a length-5 vector.")
    }
  } else {
    sds <- runif(5, y_sd_range[1], y_sd_range[2])
  }
  
  make_y <- function(name, sd_i) {
    eff <- if (name %in% dif_items) (base_loading + gamma*X) else base_loading
    eff*Lat + rnorm(N, 0, sd_i)
  }
  y1 <- make_y("y1", sds[1]); y2 <- make_y("y2", sds[2]); y3 <- make_y("y3", sds[3])
  y4 <- make_y("y4", sds[4]); y5 <- make_y("y5", sds[5])
  data.frame(Group, x1,x2,x3,x4,x5, y1,y2,y3,y4,y5, check.names = FALSE)
}


# -------- models --------
model_latX_labeled <- "
latent.model:
  latX ~ 1@0 Group;
  lat  ~ 1@0;

item.model:
  x1 ~ 1 latX@1; x2 ~ 1 latX; x3 ~ 1 latX; x4 ~ 1 latX; x5 ~ 1 latX;

  y1 ~ 1 lat@1 latX@0 latX*lat@0;
  y2 ~ 1 lat      latX     latX*lat@b_int_y2;
  y3 ~ 1 lat      latX     latX*lat@b_int_y3;
  y4 ~ 1 lat      latX     latX*lat@b_int_y4;
  y5 ~ 1 lat      latX     latX*lat@b_int_y5;
"
model_group_labeled <- "
latent.model:
  lat ~ 1@0;

item.model:
  y1 ~ 1 Group lat@1 Group*lat@0;
  y2 ~ 1 Group lat     Group*lat@b_gint_y2;
  y3 ~ 1 Group lat     Group*lat@b_gint_y3;
  y4 ~ 1 Group lat     Group*lat@b_gint_y4;
  y5 ~ 1 Group lat     Group*lat@b_gint_y5;
"
# Nulls
model_latX_null <- "
latent.model:
  latX ~ 1@0 Group;
  lat  ~ 1@0;

item.model:
  x1 ~ 1 latX@1; x2 ~ 1 latX; x3 ~ 1 latX; x4 ~ 1 latX; x5 ~ 1 latX;
  y1 ~ 1 lat@1  latX@0 latX*lat@0;
  y2 ~ 1 lat    latX@0 latX*lat@0;
  y3 ~ 1 lat    latX@0 latX*lat@0;
  y4 ~ 1 lat    latX@0 latX*lat@0;
  y5 ~ 1 lat    latX@0 latX*lat@0;
"
model_grp_null <- "
latent.model:
  lat ~ 1@0;

item.model:
  y1 ~ 1 Group@0 lat@1 Group*lat@0;
  y2 ~ 1 Group@0 lat    Group*lat@0;
  y3 ~ 1 Group@0 lat    Group*lat@0;
  y4 ~ 1 Group@0 lat    Group*lat@0;
  y5 ~ 1 Group@0 lat    Group*lat@0;
"

# -------- priors (interactions only) --------
prior_latX <- list(
  tight    = "b_int_y2 ~ normal(0, 0.25); b_int_y3 ~ normal(0, 0.25); b_int_y4 ~ normal(0, 0.25); b_int_y5 ~ normal(0, 0.25);",
  baseline = "b_int_y2 ~ normal(0, 0.50); b_int_y3 ~ normal(0, 0.50); b_int_y4 ~ normal(0, 0.50); b_int_y5 ~ normal(0, 0.50);",
  loose    = "b_int_y2 ~ normal(0, 1.00); b_int_y3 ~ normal(0, 1.00); b_int_y4 ~ normal(0, 1.00); b_int_y5 ~ normal(0, 1.00);"
)
prior_group <- list(
  tight    = "b_gint_y2 ~ normal(0, 0.25); b_gint_y3 ~ normal(0, 0.25); b_gint_y4 ~ normal(0, 0.25); b_gint_y5 ~ normal(0, 0.25);",
  baseline = "b_gint_y2 ~ normal(0, 0.50); b_gint_y3 ~ normal(0, 0.50); b_gint_y4 ~ normal(0, 0.50); b_gint_y5 ~ normal(0, 0.50);",
  loose    = "b_gint_y2 ~ normal(0, 1.00); b_gint_y3 ~ normal(0, 1.00); b_gint_y4 ~ normal(0, 1.00); b_gint_y5 ~ normal(0, 1.00);"
)

# -------- Wald specs --------
wald_tests_for <- function(model_family, dif_items) {
  prefix  <- if (model_family == "latX") "b_int_" else "b_gint_"
  singles <- paste0(prefix, c("y2","y3","y4","y5"), " = 0")
  joint_true <- paste(c(paste0(prefix, dif_items), "= 0"), collapse = " ")
  c(singles, joint_true)
}

# -------- deterministic seeds (include replication r) --------
data_seed <- function(N, dif_items, r) {
  dif_off <- if (identical(sort(dif_items), c("y2","y3"))) 10L else 20L
  as.integer(9e6 + 1000L*N + dif_off + r)
}
fit_seed <- function(family, prior, N, dif_items, r) {
  base <- switch(family,
                 "latX_full"  = 1000L,
                 "group_full" = 2000L,
                 "latX_null"  = 3000L,
                 "group_null" = 4000L,
                 5000L)
  prior_off <- if (is.null(prior)) 0L else switch(prior, tight=1L, baseline=2L, loose=3L, 0L)
  dif_off <- if (identical(sort(dif_items), c("y2","y3"))) 10L else 20L
  as.integer(base + 10L * N + dif_off + prior_off + 10000L * r)
}

# -------- small packers --------
pack_full <- function(fit, keep_fit_objects = FALSE) {
  w <- tryCatch(fit@waldtest, error = function(e) NULL)
  if (is.null(w)) w <- tryCatch(fit@tests$wald, error = function(e) NULL)
  list(
    fit            = if (keep_fit_objects) fit else NULL,
    estimates_raw  = tryCatch(fit@estimates, error = function(e) NULL),
    wald_raw       = w,
    modelfit_text  = capture.output(rblimp::modelfit(fit))
  )
}
pack_null <- function(fit, keep_fit_objects = FALSE) {
  list(
    fit            = if (keep_fit_objects) fit else NULL,
    estimates_raw  = tryCatch(fit@estimates, error = function(e) NULL),
    modelfit_text  = capture.output(rblimp::modelfit(fit))
  )
}

# -------- one replication (RAW ONLY; simplified null) --------
fit_one_rep <- function(N, dif_items, r,
                        priors_to_run, burn, iter,
                        keep_fit_objects = FALSE) {
  set.seed(data_seed(N, dif_items, r))
  dat <- gen_data(N = N, dif_items = dif_items, y_sd = Y_SD)  # <- pass fixed SD
  
  # FULL fits: Latent-X
  lat_fulls <- lapply(priors_to_run, function(pr) {
    fit <- rblimp(
      data = dat, latent = "latX lat", model = model_latX_labeled,
      burn = burn, iter = iter,
      seed = fit_seed("latX_full", pr, N, dif_items, r),
      parameters = prior_latX[[pr]],
      waldtest  = wald_tests_for("latX", dif_items),
      output    = OUTPUT_SPEC
    )
    out <- pack_full(fit, keep_fit_objects)
    if (!keep_fit_objects) rm(fit); gc(FALSE)
    out
  })
  names(lat_fulls) <- priors_to_run
  
  # FULL fits: Group
  grp_fulls <- lapply(priors_to_run, function(pr) {
    fit <- rblimp(
      data = dat, latent = "lat", model = model_group_labeled,
      burn = burn, iter = iter,
      seed = fit_seed("group_full", pr, N, dif_items, r),
      parameters = prior_group[[pr]],
      waldtest  = wald_tests_for("group", dif_items),
      output    = OUTPUT_SPEC
    )
    out <- pack_full(fit, keep_fit_objects)
    if (!keep_fit_objects) rm(fit); gc(FALSE)
    out
  })
  names(grp_fulls) <- priors_to_run
  
  # NULL fits: single canonical node per family (no wald)
  fitN <- rblimp(
    data = dat, latent = "latX lat", model = model_latX_null,
    burn = burn, iter = iter,
    seed = fit_seed("latX_null", NULL, N, dif_items, r),
    output = OUTPUT_SPEC
  )
  lat_null <- pack_null(fitN, keep_fit_objects)
  if (!keep_fit_objects) rm(fitN); gc(FALSE)
  
  fitN <- rblimp(
    data = dat, latent = "lat", model = model_grp_null,
    burn = burn, iter = iter,
    seed = fit_seed("group_null", NULL, N, dif_items, r),
    output = OUTPUT_SPEC
  )
  grp_null <- pack_null(fitN, keep_fit_objects)
  if (!keep_fit_objects) rm(fitN); gc(FALSE)
  
  # package and clean
  out_rep <- list(
    rep   = r,
    meta  = list(N = N, dif_items = dif_items, priors_run = priors_to_run),
    latX  = list(full = lat_fulls, null = lat_null),
    group = list(full = grp_fulls, null = grp_null)
  )
  rm(dat); gc(FALSE)
  out_rep
}

# -------- run simulation (each replication is [[i]]) --------
results_sim <- list()
for (N in Ns) {
  for (lbl in names(DIF_sets)) {
    dif_items <- DIF_sets[[lbl]]
    key <- sprintf("N_%04d_%s", N, lbl)
    reps <- vector("list", n_reps)
    for (r in seq_len(n_reps)) {
      reps[[r]] <- fit_one_rep(
        N, dif_items, r,
        priors_to_run, burn, iter,
        keep_fit_objects = keep_fit_objects
      )
      message("Done rep ", r, " | ", key)
    }
    results_sim[[key]] <- reps
    message("Finished scenario: ", key)
  }
}
