# ==============================================================================
# CLAIMS SEVERITY — FEATURE SELECTION COMPARISON
# Log-Normal GLM: Baseline vs Stepwise AIC vs Boruta vs Elastic Net
# ==============================================================================
#
# Compares three feature selection methods applied to the Log-Normal severity
# model. The full model (all encoded predictors) serves as the baseline.
# All methodological steps mirror the baseline exactly:
#   - Client-stratified 80/20 train-test split (seed 42)
#   - Claimant-only subsample (n_claims_year > 0)
#   - Grouped 5-fold CV (client-level, leak-free)
#   - Per-fold recipe encoding (step_novel, step_unknown, step_dummy, step_zv)
#   - Smearing back-transform: exp(mu_hat + sigma^2 / 2)
#   - Primary metric: Gini (Lorenz-curve concentration index on original cost scale)
#   - Secondary: MAE, R2_COR (squared Pearson correlation on original scale), LN_Dev
#   - RMSE tracked but not used for ranking (dominated by catastrophic tail claims)
#   - CV mean +/- SE reported throughout
#
# CV ARCHITECTURE (fold-blind selection)
# ---------------------------------------
# Each selection algorithm runs inside each CV fold on the fold's own encoded
# training data, so the validation set is never seen during selection.
#   - Stepwise AIC : full LM fit -> stepAIC inside each fold (direction = "both")
#   - Boruta       : Boruta() inside each fold on log(y_train)
#   - Elastic Net  : cv.glmnet inside each outer fold x ALPHA_GRID (5 values);
#                    best alpha selected by minimum mean outer-fold LN_Dev
# A separate full-training run of each algorithm produces the final feature
# set used for test-set evaluation and the saved production model.
#
# OUTPUTS
#   severity_feature_selection_results.csv  — raw numeric metrics (all methods)
#   severity_selected_features.csv          — selected feature sets (long format)
#   severity_feature_selection_table.html   — gt comparison table
#
# ==============================================================================


# ==============================================================================
# SECTION 0: LIBRARIES
# ==============================================================================

library(here)
library(dplyr)
library(rsample)
library(recipes)
library(MASS)
library(Boruta)
library(glmnet)
library(gt)


# ==============================================================================
# SECTION 1: DATA LOADING
# ==============================================================================

basic_glm_data       <- read.csv(here("Feature_Selection_Four_Wheels_GLM_R.csv"))
basic_glm_data_no_na <- basic_glm_data[!is.na(basic_glm_data$avg_claim_cost), ]

cat(sprintf("Full dataset: %d rows after removing NA avg_claim_cost\n",
            nrow(basic_glm_data_no_na)))


# ==============================================================================
# SECTION 2: VARIABLE DEFINITIONS
# ==============================================================================

ordinal_levels <- list(
  age_category                     = c("18-25","26-35","36-45","46-55","56-65","65+"),
  age_when_licensed_band           = c("<=25",">25"),
  driving_exp_category             = c("<=5yrs","6-10yrs","11-20yrs","20+yrs"),
  vehicle_age_category             = c("0-3yrs","4-7yrs","8-12yrs","12+yrs"),
  power_band                       = c("<100",">100"),
  cylinder_band                    = c("0-1400","1400-1600","1600-2000","2000+"),
  value_band                       = c("0-17500","17500-25000","25000-32500","32500+"),
  weight_band                      = c("0-1100","1100-1500","1500+"),
  length_imputed_band              = c("<4","4-4.5",">4.5"),
  n_doors_band                     = c("<=4",">4"),
  claims_history_band              = c("zero","one","two_to_three","four_to_five","six_plus"),
  r_claims_history_band            = c("0","0-0.35","0.35-0.7","0.7-1.5","1.5+"),
  avg_claims_per_year_history_band = c("0","0-1",">1"),
  cross_sell_ratio_band            = c("<0.75",">=0.75"),
  policy_retention_rate_band       = c("<1","1"),
  risk_score_band                  = c("0","0-0.01","0.01-0.15","0.15-0.2","0.2+"),
  premium_band                     = c("0-250","250-300","300-350","350-400","400+"),
  premium_vs_median_band           = c("-15%","-15% to 0%","0% to 15%","15% to 30%","30%+")
)

always_exclude  <- c("id","id2","n_claims_year","avg_claim_cost",
                     "lapse_flag","days_to_event","exposure_years",
                     "r_claims_history_band","claims_history_band")
premium_exclude <- c("premium_band","premium_vs_median_band")

target    <- "avg_claim_cost"
THRESHOLD <- 900


# ==============================================================================
# SECTION 3: CLIENT-STRATIFIED TRAIN-TEST SPLIT (80:20)
# ==============================================================================

set.seed(42)

client_strat <- basic_glm_data_no_na %>%
  group_by(id) %>%
  summarise(
    ever_lapsed  = as.integer(any(lapse_flag == 1,   na.rm = TRUE)),
    ever_claimed = as.integer(any(n_claims_year > 0, na.rm = TRUE)),
    n_policies   = n(),
    .groups = "drop"
  ) %>%
  mutate(
    policy_tier = case_when(
      n_policies == 1 ~ "1",
      n_policies == 2 ~ "2",
      TRUE            ~ "3+"
    ),
    strata = paste(ever_lapsed, ever_claimed, policy_tier, sep = "_")
  )

rare         <- names(table(client_strat$strata)[table(client_strat$strata) == 1])
client_strat <- client_strat %>%
  mutate(strata = if_else(strata %in% rare, "other", strata))

split_obj  <- initial_split(client_strat, prop = 0.8, strata = strata)
train_data <- basic_glm_data_no_na %>% filter(id %in% training(split_obj)$id)
test_data  <- basic_glm_data_no_na %>% filter(id %in% testing(split_obj)$id)

stopifnot(length(intersect(unique(train_data$id), unique(test_data$id))) == 0)
cat(sprintf("Train: %d policies | %d clients\nTest : %d policies | %d clients\n",
            nrow(train_data),  n_distinct(train_data$id),
            nrow(test_data),   n_distinct(test_data$id)))


# ==============================================================================
# SECTION 4: CLAIMANT SUBSAMPLES
# ==============================================================================

train_sev <- train_data %>% filter(n_claims_year > 0)
test_sev  <- test_data  %>% filter(n_claims_year > 0)

cat(sprintf("Claimant rows — Train: %d  |  Test: %d\n",
            nrow(train_sev), nrow(test_sev)))


# ==============================================================================
# SECTION 5: GROUPED 5-FOLD CV STRUCTURE
# ==============================================================================

set.seed(42)
cv_folds <- group_vfold_cv(train_data, group = id, v = 5)

stopifnot(all(vapply(cv_folds$splits, function(sp)
  length(intersect(unique(analysis(sp)$id),
                   unique(assessment(sp)$id))) == 0L, logical(1))))
cat("CV folds: 5 | Client-leakage-free: TRUE\n")


# ==============================================================================
# SECTION 6: HELPER FUNCTIONS
# ==============================================================================

recode_ordinals <- function(df) {
  for (v in intersect(names(ordinal_levels), names(df)))
    df[[v]] <- factor(df[[v]], levels = ordinal_levels[[v]], ordered = FALSE)
  df
}

get_predictor_cols <- function(data, tgt)
  setdiff(names(data), unique(c(always_exclude, premium_exclude, tgt)))

encode_predictors <- function(train_df, test_df, tgt) {
  pcols    <- get_predictor_cols(train_df, tgt)
  x_tr     <- recode_ordinals(train_df[, pcols, drop = FALSE])
  x_te     <- recode_ordinals(test_df[,  pcols, drop = FALSE])
  rec      <- recipe(~., data = x_tr) %>%
    step_novel(all_nominal_predictors())   %>%
    step_unknown(all_nominal_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_zv(all_predictors())
  prep_rec <- prep(rec, training = x_tr, retain = TRUE)
  list(
    train  = bake(prep_rec, new_data = NULL),
    test   = bake(prep_rec, new_data = x_te),
    recipe = prep_rec
  )
}

# Embeds the data frame directly into the lm() call via bquote/eval so that
# stepAIC's internal update() can re-evaluate the call without needing to
# resolve 'y' or 'X' from a long-gone execution frame.
fit_lm <- function(y, X) {
  df <- data.frame(resp__ = y, X)
  eval(bquote(lm(resp__ ~ ., data = .(df))))
}

rmse       <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae        <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
r2_cor     <- function(obs, pred) cor(obs, pred, use = "complete.obs")^2

gini_index <- function(obs, pred) {
  n   <- length(obs)
  ord <- order(pred, decreasing = TRUE)
  co  <- c(0, cumsum(obs[ord]) / sum(obs))
  cp  <- seq(0, 1, length.out = n + 1)
  2 * sum(diff(cp) * (co[-1] + co[-(n + 1)]) / 2) - 1
}

ln_backtr    <- function(pred_log, sigma) exp(pred_log + sigma^2 / 2)

ll_lognormal <- function(y, mu_log, sigma)
  sum(dnorm(log(y), mean = mu_log, sd = sigma, log = TRUE) - log(y))

ln_aic_bic <- function(y, mu_log, sigma, p, n) {
  ll <- ll_lognormal(y, mu_log, sigma)
  k  <- p + 1
  c(AIC = -2 * ll + 2 * k, BIC = -2 * ll + log(n) * k)
}

sev_metrics <- function(y, pred, r2, aic, bic, pred_log = NULL) {
  ln_dev_val <- if (!is.null(pred_log))
    sum((log(y) - pred_log)^2, na.rm = TRUE)
  else
    sum((log(y) - log(pred))^2, na.rm = TRUE)
  c(R2     = round(r2, 4),
    RMSE   = round(rmse(y, pred), 2),
    MAE    = round(mae(y, pred),  2),
    LN_Dev = round(ln_dev_val,    2),
    Gini   = round(gini_index(y, pred), 4),
    AIC    = round(aic, 0),
    BIC    = round(bic, 0))
}

cv_summary <- function(tbl) {
  if (nrow(tbl) == 0L) {
    na_vec <- setNames(rep(NA_real_, ncol(tbl)), names(tbl))
    return(list(mean = na_vec, se = na_vec))
  }
  list(
    mean = colMeans(tbl, na.rm = TRUE),
    se   = apply(tbl, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  )
}

fmt_cv <- function(m, s, digits = 4)
  sprintf(paste0("%.", digits, "f (\u00b1%.", digits, "f)"), m, s)


# ==============================================================================
# SECTION 7: FULL-CLAIMANT ENCODING (SHARED ACROSS ALL MODELS)
# ==============================================================================

enc_sev  <- encode_predictors(train_sev, test_sev, target)
X_tr_sev <- enc_sev$train
X_te_sev <- enc_sev$test
y_tr_sev <- train_sev[[target]]
y_te_sev <- test_sev[[target]]
n_tr_sev <- length(y_tr_sev)
n_te_sev <- length(y_te_sev)

cat(sprintf("Encoded predictor columns: %d\n", ncol(X_tr_sev)))


# ==============================================================================
# SECTION 8: FOLD ID VECTOR FOR ELASTIC NET
# ==============================================================================

fold_id_vec <- integer(nrow(train_sev))
for (k in seq_along(cv_folds$splits)) {
  assess_ids               <- assessment(cv_folds$splits[[k]])$id
  fold_id_vec[train_sev$id %in% assess_ids] <- k
}

stopifnot(all(fold_id_vec > 0))
cat(sprintf("Fold ID vector: %s (fold sizes)\n",
            paste(table(fold_id_vec), collapse = " | ")))


# ==============================================================================
# SECTION 9: REUSABLE CV AND TEST EVALUATION FUNCTIONS
# ==============================================================================

# run_cv_ln: used by the baseline only. For the three selection methods,
# selection runs inside each fold (see per-fold lapply loops in Sections 11-13).

run_cv_ln <- function(feature_set,
                      cv_folds,
                      train_data,
                      target  = "avg_claim_cost",
                      verbose = TRUE) {
  n_folds <- length(cv_folds$splits)

  results <- lapply(seq_along(cv_folds$splits), function(i) {
    sp     <- cv_folds$splits[[i]]
    t_fold <- proc.time()[["elapsed"]]
    if (verbose) cat(sprintf("    Fold %d/%d ... ", i, n_folds))

    tr <- analysis(sp)   %>% filter(n_claims_year > 0)
    va <- assessment(sp) %>% filter(n_claims_year > 0)
    if (nrow(tr) == 0 || nrow(va) == 0) {
      if (verbose) cat("skipped (no claimant rows)\n")
      return(NULL)
    }

    enc <- encode_predictors(tr, va, target)
    ytr <- tr[[target]]
    yva <- va[[target]]

    sel <- intersect(feature_set,
                     intersect(names(enc$train), names(enc$test)))
    if (length(sel) == 0) {
      if (verbose) cat("skipped (no features survive ZV filter)\n")
      return(NULL)
    }

    m <- tryCatch(
      fit_lm(log(ytr), enc$train[, sel, drop = FALSE]),
      error = function(e) NULL
    )
    if (is.null(m)) {
      if (verbose) cat("skipped (lm failed)\n")
      return(NULL)
    }

    sg       <- sigma(m)
    pred_log <- predict(m, enc$test[, sel, drop = FALSE])
    pred     <- ln_backtr(pred_log, sg)
    k        <- length(coef(m))
    n        <- length(yva)
    ab       <- ln_aic_bic(yva, pred_log, sg, k, n)

    elapsed <- proc.time()[["elapsed"]] - t_fold
    if (verbose) cat(sprintf("done (%.1fs)  Gini=%.4f  MAE=%.1f\n",
                             elapsed, gini_index(yva, pred), mae(yva, pred)))

    data.frame(
      R2     = r2_cor(yva, pred),
      RMSE   = rmse(yva, pred),
      MAE    = mae(yva, pred),
      LN_Dev = sum((log(yva) - pred_log)^2),
      Gini   = gini_index(yva, pred),
      AIC    = ab[["AIC"]],
      BIC    = ab[["BIC"]]
    )
  })

  bind_rows(Filter(Negate(is.null), results))
}


run_test_ln <- function(feature_set, X_tr, y_tr, X_te, y_te) {
  sel <- intersect(feature_set,
                   intersect(names(X_tr), names(X_te)))

  if (length(sel) == 0L)
    stop("run_test_ln: feature_set is empty after intersection with encoded ",
         "column names. The upstream selection step produced no usable features.")

  fit      <- fit_lm(log(y_tr), X_tr[, sel, drop = FALSE])
  sg       <- sigma(fit)
  pred_log <- predict(fit, X_te[, sel, drop = FALSE])
  pred     <- ln_backtr(pred_log, sg)
  p        <- length(coef(fit))
  n        <- length(y_te)
  ab       <- ln_aic_bic(y_te, pred_log, sg, p, n)

  list(
    fit        = fit,
    metrics    = sev_metrics(y_te, pred, r2_cor(y_te, pred),
                             ab[["AIC"]], ab[["BIC"]], pred_log),
    n_features = length(sel)
  )
}


build_result_row <- function(method, n_features, cv_df, test_metrics) {
  s <- cv_summary(cv_df)
  data.frame(
    Method     = method,
    N_Features = n_features,
    CV_Gini    = fmt_cv(s$mean["Gini"],   s$se["Gini"],   4),
    CV_MAE     = fmt_cv(s$mean["MAE"],    s$se["MAE"],    1),
    CV_R2      = fmt_cv(s$mean["R2"],     s$se["R2"],     4),
    CV_LN_Dev  = fmt_cv(s$mean["LN_Dev"], s$se["LN_Dev"], 0),
    CV_RMSE    = fmt_cv(s$mean["RMSE"],   s$se["RMSE"],   1),
    Test_Gini  = as.numeric(round(unname(test_metrics["Gini"]), 4)),
    Test_MAE   = as.numeric(round(unname(test_metrics["MAE"]),  1)),
    Test_R2    = as.numeric(round(unname(test_metrics["R2"]),   4)),
    Test_RMSE  = as.numeric(round(unname(test_metrics["RMSE"]), 1)),
    stringsAsFactors = FALSE
  )
}


make_raw_row <- function(method, n_features, cv_df, test_metrics) {
  s <- cv_summary(cv_df)
  data.frame(
    Method        = method,
    N_Features    = n_features,
    CV_Gini_mean  = round(unname(s$mean["Gini"]),   4),
    CV_Gini_se    = round(unname(s$se["Gini"]),     4),
    CV_MAE_mean   = round(unname(s$mean["MAE"]),    2),
    CV_MAE_se     = round(unname(s$se["MAE"]),      2),
    CV_R2_mean    = round(unname(s$mean["R2"]),     4),
    CV_R2_se      = round(unname(s$se["R2"]),       4),
    CV_LNDev_mean = round(unname(s$mean["LN_Dev"]), 0),
    CV_LNDev_se   = round(unname(s$se["LN_Dev"]),   0),
    CV_RMSE_mean  = round(unname(s$mean["RMSE"]),   2),
    CV_RMSE_se    = round(unname(s$se["RMSE"]),     2),
    Test_Gini     = round(unname(test_metrics["Gini"]),   4),
    Test_MAE      = round(unname(test_metrics["MAE"]),    2),
    Test_R2       = round(unname(test_metrics["R2"]),     4),
    Test_LNDev    = round(unname(test_metrics["LN_Dev"]), 0),
    Test_RMSE     = round(unname(test_metrics["RMSE"]),   2),
    stringsAsFactors = FALSE
  )
}


save_method_results <- function(file_stem, method_label,
                                sel_features, cv_raw, test_result) {
  bundle <- list(
    selected_features = sel_features,
    cv_raw            = cv_raw,
    cv_summary        = cv_summary(cv_raw),
    test_metrics      = test_result$metrics,
    model_fit         = test_result$fit,
    n_features        = test_result$n_features
  )
  rds_path <- here(sprintf("severity_fs_%s.rds",         file_stem))
  csv_path <- here(sprintf("severity_fs_%s_metrics.csv", file_stem))
  saveRDS(bundle, rds_path)
  write.csv(
    make_raw_row(method_label, test_result$n_features, cv_raw, test_result$metrics),
    csv_path, row.names = FALSE
  )
  cat(sprintf("  [Saved] %-40s | %s\n", basename(rds_path), basename(csv_path)))
  invisible(bundle)
}


# ==============================================================================
# FEATURE SELECTION CONSTANTS
# ==============================================================================

BORUTA_MAXRUNS <- 150
ALPHA_GRID     <- c(0.1, 0.25, 0.5, 0.75, 1.0)


# ==============================================================================
# SECTION 10: BASELINE — FULL LOG-NORMAL MODEL
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("BASELINE: FULL LOG-NORMAL MODEL\n")
cat(strrep("=", 70), "\n")

baseline_features <- names(X_tr_sev)
cat(sprintf("Baseline uses all %d encoded predictors\n", length(baseline_features)))

t_start <- Sys.time()
cat("Running 5-fold CV...\n")
cv_baseline <- run_cv_ln(baseline_features, cv_folds, train_data)
cv_s_base   <- cv_summary(cv_baseline)
cat(sprintf("CV complete (%.1f mins)\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

cat(sprintf("  CV Gini:  %.4f (+/- %.4f)\n", cv_s_base$mean["Gini"],   cv_s_base$se["Gini"]))
cat(sprintf("  CV MAE:   %.1f  (+/- %.1f)\n",  cv_s_base$mean["MAE"],    cv_s_base$se["MAE"]))
cat(sprintf("  CV R2:    %.4f (+/- %.4f)\n", cv_s_base$mean["R2"],     cv_s_base$se["R2"]))
cat(sprintf("  CV LNDev: %.0f  (+/- %.0f)\n",  cv_s_base$mean["LN_Dev"], cv_s_base$se["LN_Dev"]))

cat("Evaluating on hold-out test set...\n")
test_baseline <- run_test_ln(baseline_features, X_tr_sev, y_tr_sev, X_te_sev, y_te_sev)
cat("Test metrics (Baseline):\n")
print(test_baseline$metrics)

save_method_results("baseline", "Baseline (Full)",
                    baseline_features, cv_baseline, test_baseline)

row_baseline     <- build_result_row("Baseline (Full)", test_baseline$n_features,
                                     cv_baseline, test_baseline$metrics)
raw_row_baseline <- make_raw_row("Baseline (Full)", test_baseline$n_features,
                                 cv_baseline, test_baseline$metrics)


# ==============================================================================
# SECTION 11: FEATURE SELECTION — STEPWISE AIC
# ==============================================================================
#
# stepAIC runs inside each CV fold on that fold's own encoded training data
# (direction = "both", trace = FALSE). The fold's full LM provides the
# starting point. Selected features are the non-intercept coefficient names
# from the stepAIC output; the LM is then refit on those features to produce
# fold-level metrics.
#
# Full training run: stepAIC on all training claimants (ln_fit_full as start)
# -> selected_step -> test evaluation and saved model.

cat("\n", strrep("=", 70), "\n")
cat("FEATURE SELECTION: STEPWISE AIC\n")
cat(strrep("=", 70), "\n")

t_start <- Sys.time()
cat("Running 5-fold CV (stepAIC inside each fold)...\n")

cv_step_raw <- lapply(seq_along(cv_folds$splits), function(i) {
  sp     <- cv_folds$splits[[i]]
  t_fold <- proc.time()[["elapsed"]]
  cat(sprintf("    Fold %d/%d ... ", i, length(cv_folds$splits)))

  tr <- analysis(sp)   %>% filter(n_claims_year > 0)
  va <- assessment(sp) %>% filter(n_claims_year > 0)
  if (nrow(tr) == 0 || nrow(va) == 0) {
    cat("skipped (no claimant rows)\n"); return(NULL)
  }

  enc <- encode_predictors(tr, va, target)
  ytr <- tr[[target]]; yva <- va[[target]]

  full_m <- tryCatch(
    fit_lm(log(ytr), enc$train),
    error = function(e) NULL
  )
  if (is.null(full_m)) { cat("skipped (full LM failed)\n"); return(NULL) }

  step_m <- tryCatch(
    MASS::stepAIC(full_m, direction = "both", trace = FALSE),
    error = function(e) NULL
  )
  if (is.null(step_m)) { cat("skipped (stepAIC failed)\n"); return(NULL) }

  sel <- intersect(
    setdiff(names(coef(step_m)), "(Intercept)"),
    intersect(names(enc$train), names(enc$test))
  )
  if (length(sel) == 0) { cat("skipped (no features survived)\n"); return(NULL) }

  m <- tryCatch(
    fit_lm(log(ytr), enc$train[, sel, drop = FALSE]),
    error = function(e) NULL
  )
  if (is.null(m)) { cat("skipped (refit failed)\n"); return(NULL) }

  sg       <- sigma(m)
  pred_log <- predict(m, enc$test[, sel, drop = FALSE])
  pred     <- ln_backtr(pred_log, sg)
  k        <- length(coef(m)); n <- length(yva)
  ab       <- ln_aic_bic(yva, pred_log, sg, k, n)

  elapsed <- proc.time()[["elapsed"]] - t_fold
  cat(sprintf("done (%.1fs)  vars=%d  Gini=%.4f  MAE=%.1f\n",
              elapsed, length(sel), gini_index(yva, pred), mae(yva, pred)))

  data.frame(
    R2     = r2_cor(yva, pred),
    RMSE   = rmse(yva, pred),
    MAE    = mae(yva, pred),
    LN_Dev = sum((log(yva) - pred_log)^2),
    Gini   = gini_index(yva, pred),
    AIC    = ab[["AIC"]],
    BIC    = ab[["BIC"]]
  )
})

cv_step   <- bind_rows(Filter(Negate(is.null), cv_step_raw))
cv_s_step <- cv_summary(cv_step)
cat(sprintf("CV complete (%.1f mins)\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat(sprintf("  CV Gini:  %.4f (+/- %.4f)\n", cv_s_step$mean["Gini"],   cv_s_step$se["Gini"]))
cat(sprintf("  CV MAE:   %.1f  (+/- %.1f)\n",  cv_s_step$mean["MAE"],    cv_s_step$se["MAE"]))
cat(sprintf("  CV R2:    %.4f (+/- %.4f)\n", cv_s_step$mean["R2"],     cv_s_step$se["R2"]))
cat(sprintf("  CV LNDev: %.0f  (+/- %.0f)\n",  cv_s_step$mean["LN_Dev"], cv_s_step$se["LN_Dev"]))

ln_fit_full <- fit_lm(log(y_tr_sev), X_tr_sev)
t_start     <- Sys.time()
cat("Running stepAIC on full training set (direction = 'both')...\n")
ln_step       <- MASS::stepAIC(ln_fit_full, direction = "both", trace = FALSE)
selected_step <- setdiff(names(coef(ln_step)), "(Intercept)")
cat(sprintf("stepAIC complete (%.1f mins) — selected %d features from %d\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins")),
            length(selected_step), ncol(X_tr_sev)))

cat("Evaluating on hold-out test set...\n")
test_step <- run_test_ln(selected_step, X_tr_sev, y_tr_sev, X_te_sev, y_te_sev)
cat("Test metrics (Stepwise AIC):\n"); print(test_step$metrics)

save_method_results("stepwise", "Stepwise AIC", selected_step, cv_step, test_step)

row_step     <- build_result_row("Stepwise AIC", length(selected_step),
                                 cv_step, test_step$metrics)
raw_row_step <- make_raw_row("Stepwise AIC", length(selected_step),
                             cv_step, test_step$metrics)


# ==============================================================================
# SECTION 12: FEATURE SELECTION — BORUTA
# ==============================================================================
#
# Boruta runs inside each CV fold on that fold's own encoded log(y) training
# data (fold-blind selection). Hyperparameters are consistent across all calls:
#   num.trees    = 500           : sufficient for importance stability on ~15K rows
#   maxRuns      = BORUTA_MAXRUNS: 150, consistent with frequency pipeline
#   doTrace      = 2             : numbered per-run progress lines
#   verbose      = FALSE         : suppresses ranger's per-tree output
#
# After Boruta confirms features, TentativeRoughFix resolves any tentatives by
# comparing their median importance against the shadow maximum. Skipping folds
# where Boruta confirms nothing is preferable to falling back to all features,
# which would make those folds indistinguishable from the baseline and bias the
# CV mean upward.

cat("\n", strrep("=", 70), "\n")
cat("FEATURE SELECTION: BORUTA\n")
cat(strrep("=", 70), "\n")

t_start <- Sys.time()
cat(sprintf("Running 5-fold CV (Boruta inside each fold, maxRuns=%d)...\n",
            BORUTA_MAXRUNS))

cv_boruta_raw <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  cat(sprintf("\n  \u2500\u2500 Fold %d/%d \u2014 running Boruta \u2500\u2500\n",
              i, length(cv_folds$splits)))
  t_fold <- proc.time()[["elapsed"]]

  tr <- analysis(sp)   %>% filter(n_claims_year > 0)
  va <- assessment(sp) %>% filter(n_claims_year > 0)
  if (nrow(tr) == 0 || nrow(va) == 0) {
    cat("  skipped (no claimant rows)\n"); return(NULL)
  }

  enc <- encode_predictors(tr, va, target)
  ytr <- tr[[target]]; yva <- va[[target]]

  bor <- Boruta::Boruta(
    x         = enc$train,
    y         = log(ytr),
    num.trees = 500,
    maxRuns   = BORUTA_MAXRUNS,
    doTrace   = 2,
    verbose   = FALSE
  )
  bor       <- Boruta::TentativeRoughFix(bor)
  confirmed <- Boruta::getSelectedAttributes(bor, withTentative = FALSE)
  if (length(confirmed) == 0L) {
    cat("  skipped (Boruta confirmed no features after TentativeRoughFix)\n")
    return(NULL)
  }

  sel <- intersect(confirmed, intersect(names(enc$train), names(enc$test)))
  if (length(sel) == 0) { cat("  skipped (no features survived)\n"); return(NULL) }

  m <- tryCatch(
    fit_lm(log(ytr), enc$train[, sel, drop = FALSE]),
    error = function(e) NULL
  )
  if (is.null(m)) { cat("  skipped (LM failed)\n"); return(NULL) }

  sg       <- sigma(m)
  pred_log <- predict(m, enc$test[, sel, drop = FALSE])
  pred     <- ln_backtr(pred_log, sg)
  k        <- length(coef(m)); n <- length(yva)
  ab       <- ln_aic_bic(yva, pred_log, sg, k, n)

  elapsed <- proc.time()[["elapsed"]] - t_fold
  cat(sprintf("\n  -> Fold %d/%d done (%.1f mins) \u2014 %d confirmed | Gini=%.4f  MAE=%.1f\n",
              i, length(cv_folds$splits), elapsed / 60,
              length(sel), gini_index(yva, pred), mae(yva, pred)))

  data.frame(
    R2     = r2_cor(yva, pred),
    RMSE   = rmse(yva, pred),
    MAE    = mae(yva, pred),
    LN_Dev = sum((log(yva) - pred_log)^2),
    Gini   = gini_index(yva, pred),
    AIC    = ab[["AIC"]],
    BIC    = ab[["BIC"]]
  )
})

cv_boruta   <- bind_rows(Filter(Negate(is.null), cv_boruta_raw))
cv_s_boruta <- cv_summary(cv_boruta)
cat(sprintf("\nCV complete (%.1f mins)\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat(sprintf("  CV Gini:  %.4f (+/- %.4f)\n", cv_s_boruta$mean["Gini"],   cv_s_boruta$se["Gini"]))
cat(sprintf("  CV MAE:   %.1f  (+/- %.1f)\n",  cv_s_boruta$mean["MAE"],    cv_s_boruta$se["MAE"]))
cat(sprintf("  CV R2:    %.4f (+/- %.4f)\n", cv_s_boruta$mean["R2"],     cv_s_boruta$se["R2"]))
cat(sprintf("  CV LNDev: %.0f  (+/- %.0f)\n",  cv_s_boruta$mean["LN_Dev"], cv_s_boruta$se["LN_Dev"]))

set.seed(42)
t_start <- Sys.time()
cat(sprintf("\nRunning Boruta on full training set (num.trees=500, maxRuns=%d)...\n",
            BORUTA_MAXRUNS))

boruta_fit <- Boruta::Boruta(
  x         = X_tr_sev,
  y         = log(y_tr_sev),
  num.trees = 500,
  maxRuns   = BORUTA_MAXRUNS,
  doTrace   = 2,
  verbose   = FALSE
)

cat(sprintf("\nBoruta complete (%.1f mins)\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

boruta_fixed    <- Boruta::TentativeRoughFix(boruta_fit)
selected_boruta <- Boruta::getSelectedAttributes(boruta_fixed)
if (length(selected_boruta) == 0L)
  stop("Boruta rejected all features on the full training set. ",
       "Consider increasing maxRuns or inspecting the importance distribution.")
cat(sprintf("Selected %d features from %d (removed %d)\n",
            length(selected_boruta), ncol(X_tr_sev),
            ncol(X_tr_sev) - length(selected_boruta)))
cat(sprintf("  Confirmed: %d | Tentative (resolved): %d | Rejected: %d\n",
            sum(boruta_fit$finalDecision == "Confirmed"),
            sum(boruta_fit$finalDecision == "Tentative"),
            sum(boruta_fit$finalDecision == "Rejected")))

cat("Evaluating on hold-out test set...\n")
test_boruta <- run_test_ln(selected_boruta, X_tr_sev, y_tr_sev, X_te_sev, y_te_sev)
cat("Test metrics (Boruta):\n"); print(test_boruta$metrics)

save_method_results("boruta", "Boruta", selected_boruta, cv_boruta, test_boruta)

row_boruta     <- build_result_row("Boruta", length(selected_boruta),
                                   cv_boruta, test_boruta$metrics)
raw_row_boruta <- make_raw_row("Boruta", length(selected_boruta),
                               cv_boruta, test_boruta$metrics)


# ==============================================================================
# SECTION 13: FEATURE SELECTION — ELASTIC NET
# ==============================================================================
#
# Outer CV loop x ALPHA_GRID (fold-blind selection):
#   For each fold x alpha: cv.glmnet (Gaussian, inner 5-fold, lambda.1se) on
#   that fold's training claimants; predictions on the fold's validation set.
#   Sigma for back-transform estimated using a df-corrected estimator (n - p - 1)
#   to match sigma() from lm, keeping the smearing back-transform comparable
#   across all four methods.
#   Best alpha = minimum mean outer-fold LN_Dev.
#
# Full training run: cv.glmnet with best_alpha and fold_id_vec (client-grouped
#   folds for lambda selection). Non-zero coefficients at lambda.1se define
#   selected_enet. Unpenalised LM refit for test evaluation so that comparisons
#   with Stepwise AIC and Boruta remain on the same footing.

cat("\n", strrep("=", 70), "\n")
cat("FEATURE SELECTION: ELASTIC NET\n")
cat(strrep("=", 70), "\n")

t_start <- Sys.time()
cat(sprintf("Running 5-fold CV over alpha grid {%s}...\n",
            paste(ALPHA_GRID, collapse = ", ")))

cv_enet_all <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  cat(sprintf("    Fold %d/%d ... ", i, length(cv_folds$splits)))

  tr <- analysis(sp)   %>% filter(n_claims_year > 0)
  va <- assessment(sp) %>% filter(n_claims_year > 0)
  if (nrow(tr) == 0 || nrow(va) == 0) { cat("skipped\n"); return(NULL) }

  enc    <- encode_predictors(tr, va, target)
  ytr    <- tr[[target]]; yva <- va[[target]]
  X_tr_m <- as.matrix(enc$train)
  X_va_m <- as.matrix(enc$test)

  fold_results <- lapply(ALPHA_GRID, function(a) {
    set.seed(42)
    cv_fit <- tryCatch(
      glmnet::cv.glmnet(
        x           = X_tr_m,
        y           = log(ytr),
        family      = "gaussian",
        alpha       = a,
        standardize = TRUE,
        nfolds      = 5
      ),
      error = function(e) NULL
    )
    if (is.null(cv_fit)) return(NULL)

    pred_log  <- as.vector(predict(cv_fit, newx = X_va_m, s = "lambda.1se"))
    train_log <- as.vector(predict(cv_fit, newx = X_tr_m, s = "lambda.1se"))
    n_nonzero <- sum(as.vector(coef(cv_fit, s = "lambda.1se"))[-1] != 0)
    df_resid  <- max(length(ytr) - n_nonzero - 1L, 1L)
    sg        <- sqrt(sum((log(ytr) - train_log)^2) / df_resid)
    pred      <- ln_backtr(pred_log, sg)

    data.frame(
      alpha  = a,
      R2     = r2_cor(yva, pred),
      RMSE   = rmse(yva, pred),
      MAE    = mae(yva, pred),
      LN_Dev = sum((log(yva) - pred_log)^2),
      Gini   = gini_index(yva, pred),
      n_vars = n_nonzero
    )
  })

  res <- bind_rows(Filter(Negate(is.null), fold_results))
  if (nrow(res) == 0) { cat("skipped\n"); return(NULL) }
  best_row <- res[which.min(res$LN_Dev), ]
  cat(sprintf("done  best alpha=%.2f  Gini=%.4f  MAE=%.1f  vars=%d\n",
              best_row$alpha, best_row$Gini, best_row$MAE,
              as.integer(best_row$n_vars)))
  res
})

cv_enet_df <- bind_rows(Filter(Negate(is.null), cv_enet_all))
cat(sprintf("CV complete (%.1f mins)\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

enet_alpha_summary <- cv_enet_df %>%
  group_by(alpha) %>%
  summarise(
    LNDev_mean = mean(LN_Dev, na.rm = TRUE),
    LNDev_se   = sd(LN_Dev,   na.rm = TRUE) / sqrt(n()),
    Gini_mean  = mean(Gini,   na.rm = TRUE),
    MAE_mean   = mean(MAE,    na.rm = TRUE),
    nvar_mean  = mean(n_vars, na.rm = TRUE),
    .groups    = "drop"
  )
cat("\nElastic Net — CV summary by alpha:\n")
print(as.data.frame(enet_alpha_summary), digits = 4)

best_alpha <- enet_alpha_summary$alpha[which.min(enet_alpha_summary$LNDev_mean)]
cat(sprintf("\nSelected alpha: %.2f (minimum mean outer-fold LN_Dev)\n", best_alpha))

cv_enet   <- cv_enet_df %>%
  filter(alpha == best_alpha) %>%
  dplyr::select(R2, RMSE, MAE, LN_Dev, Gini)
cv_s_enet <- cv_summary(cv_enet)

cat(sprintf("  CV Gini:  %.4f (+/- %.4f)\n", cv_s_enet$mean["Gini"],   cv_s_enet$se["Gini"]))
cat(sprintf("  CV MAE:   %.1f  (+/- %.1f)\n",  cv_s_enet$mean["MAE"],    cv_s_enet$se["MAE"]))
cat(sprintf("  CV R2:    %.4f (+/- %.4f)\n", cv_s_enet$mean["R2"],     cv_s_enet$se["R2"]))
cat(sprintf("  CV LNDev: %.0f  (+/- %.0f)\n",  cv_s_enet$mean["LN_Dev"], cv_s_enet$se["LN_Dev"]))

set.seed(42)
t_start <- Sys.time()
cat(sprintf("\nRunning cv.glmnet on full training set (alpha=%.2f, lambda.1se)...\n",
            best_alpha))

enet_cv <- glmnet::cv.glmnet(
  x           = as.matrix(X_tr_sev),
  y           = log(y_tr_sev),
  alpha       = best_alpha,
  foldid      = fold_id_vec,
  standardize = TRUE,
  family      = "gaussian"
)

cat(sprintf("cv.glmnet complete (%.1f mins) \u2014 lambda.1se=%.5f\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins")),
            enet_cv$lambda.1se))

coef_1se      <- coef(enet_cv, s = "lambda.1se")
nonzero_idx   <- which(coef_1se[-1, 1] != 0)
selected_enet <- colnames(X_tr_sev)[nonzero_idx]
cat(sprintf("Selected %d features from %d (removed %d)\n",
            length(selected_enet), ncol(X_tr_sev),
            ncol(X_tr_sev) - length(selected_enet)))

cat("Evaluating on hold-out test set...\n")
test_enet <- run_test_ln(selected_enet, X_tr_sev, y_tr_sev, X_te_sev, y_te_sev)
cat("Test metrics (Elastic Net):\n"); print(test_enet$metrics)

enet_label <- sprintf("Elastic Net (a=%.2f)", best_alpha)
save_method_results("enet", enet_label, selected_enet, cv_enet, test_enet)

row_enet     <- build_result_row(enet_label, length(selected_enet),
                                 cv_enet, test_enet$metrics)
raw_row_enet <- make_raw_row(enet_label, length(selected_enet),
                             cv_enet, test_enet$metrics)


# ==============================================================================
# SECTION 14: COMPILE RESULTS AND EXPORT COMBINED OUTPUTS
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("RESULTS COMPILATION\n")
cat(strrep("=", 70), "\n")

results_table <- bind_rows(row_baseline, row_step, row_boruta, row_enet)
cat("\nFormatted summary table:\n")
print(results_table)

raw_results <- bind_rows(raw_row_baseline, raw_row_step, raw_row_boruta, raw_row_enet)

write.csv(raw_results,
          here("severity_feature_selection_results.csv"),
          row.names = FALSE)
cat("\nRaw metrics exported to: severity_feature_selection_results.csv\n")

selected_features_df <- data.frame(
  Feature = c(baseline_features, selected_step, selected_boruta, selected_enet),
  Method  = c(
    rep("Baseline",     length(baseline_features)),
    rep("Stepwise_AIC", length(selected_step)),
    rep("Boruta",       length(selected_boruta)),
    rep("Elastic_Net",  length(selected_enet))
  ),
  stringsAsFactors = FALSE
)

write.csv(selected_features_df,
          here("severity_selected_features.csv"),
          row.names = FALSE)
cat("Selected feature sets exported to: severity_selected_features.csv\n")


# ==============================================================================
# SECTION 15: FINAL COMPARISON TABLE (gt)
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("GENERATING gt COMPARISON TABLE\n")
cat(strrep("=", 70), "\n")

display_tbl <- results_table %>%
  dplyr::select(Method, N_Features,
                CV_Gini, CV_MAE, CV_R2, CV_LN_Dev,
                Test_Gini, Test_MAE, Test_R2)

best_test_gini <- max(display_tbl$Test_Gini)
best_test_mae  <- min(display_tbl$Test_MAE)
best_test_r2   <- max(display_tbl$Test_R2)

gt_table <- display_tbl %>%
  gt() %>%

  tab_header(
    title    = md("**Claims Severity: Log-Normal Feature Selection Comparison**"),
    subtitle = md(paste0(
      "*5-fold grouped client-level CV (mean \u00b1 SE) and hold-out test set",
      " &nbsp;|&nbsp; Primary metric: Gini*"))
  ) %>%

  tab_spanner(
    label   = md("**Cross-Validation (mean \u00b1 SE)**"),
    columns = c(CV_Gini, CV_MAE, CV_R2, CV_LN_Dev)
  ) %>%
  tab_spanner(
    label   = md("**Hold-out Test**"),
    columns = c(Test_Gini, Test_MAE, Test_R2)
  ) %>%

  cols_label(
    Method     = md("**Model**"),
    N_Features = md("**N Features**"),
    CV_Gini    = md("**Gini \u2191**"),
    CV_MAE     = md("**MAE (\u00a3) \u2193**"),
    CV_R2      = md("**R\u00b2 \u2191**"),
    CV_LN_Dev  = md("**LN Dev \u2193**"),
    Test_Gini  = md("**Gini \u2191**"),
    Test_MAE   = md("**MAE (\u00a3) \u2193**"),
    Test_R2    = md("**R\u00b2 \u2191**")
  ) %>%

  tab_style(
    style     = cell_fill(color = "#EFEFEF"),
    locations = cells_body(rows = Method == "Baseline (Full)")
  ) %>%

  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = Method)
  ) %>%

  tab_style(
    style     = cell_fill(color = "#C6EFCE"),
    locations = cells_body(columns = Test_Gini, rows = Test_Gini == best_test_gini)
  ) %>%

  tab_style(
    style     = cell_fill(color = "#C6EFCE"),
    locations = cells_body(columns = Test_MAE, rows = Test_MAE == best_test_mae)
  ) %>%

  tab_style(
    style     = cell_fill(color = "#C6EFCE"),
    locations = cells_body(columns = Test_R2, rows = Test_R2 == best_test_r2)
  ) %>%

  tab_footnote(
    footnote  = paste0(
      "Gini: Lorenz-curve concentration index computed on the original cost scale. ",
      "Higher = better claimant ranking. Primary metric for model selection."),
    locations = cells_column_labels(columns = CV_Gini)
  ) %>%
  tab_footnote(
    footnote  = paste0(
      "R\u00b2 = cor(y, \u0177)\u00b2 on the original cost scale for all models ",
      "and splits. Not the OLS R\u00b2 on the log scale."),
    locations = cells_column_labels(columns = CV_R2)
  ) %>%
  tab_footnote(
    footnote  = paste0(
      "LN Dev = \u03a3(log y \u2212 log \u0177\u209c\u1d63\u1d9c)\u00b2 using raw ",
      "log-scale predictions before back-transformation. ",
      "The quantity the LM directly minimised."),
    locations = cells_column_labels(columns = CV_LN_Dev)
  ) %>%
  tab_footnote(
    footnote  = paste0(
      "RMSE excluded from primary comparison: dominated by catastrophic tail claims ",
      "that no severity model recovers well."),
    locations = cells_title(groups = "subtitle")
  ) %>%

  tab_options(
    table.font.size              = px(13),
    heading.title.font.size      = px(15),
    heading.subtitle.font.size   = px(12),
    column_labels.font.weight    = "bold",
    table.border.top.color       = "#333333",
    table.border.bottom.color    = "#333333",
    heading.border.bottom.color  = "#333333",
    column_labels.border.top.color    = "#333333",
    column_labels.border.bottom.color = "#333333",
    row.striping.include_table_body   = FALSE
  )

print(gt_table)

gt::gtsave(gt_table, here("severity_feature_selection_table.html"))
cat("Table saved to: severity_feature_selection_table.html\n")

cat("\n", strrep("=", 70), "\n")
cat("Script complete.\n")
cat(strrep("=", 70), "\n")


# ==============================================================================
# SECTION 16: EXPORT COEFFICIENT TABLES FROM SAVED BUNDLES
# ==============================================================================

library(here)

baseline <- readRDS(here("severity_fs_baseline.rds"))
stepwise <- readRDS(here("severity_fs_stepwise.rds"))
boruta   <- readRDS(here("severity_fs_boruta.rds"))
enet     <- readRDS(here("severity_fs_enet.rds"))

export_glm_output <- function(model_obj, file_name) {
  coef_table          <- as.data.frame(summary(model_obj)$coefficients)
  coef_table$Variable <- rownames(coef_table)
  rownames(coef_table) <- NULL
  coef_table <- coef_table[, c("Variable", "Estimate", "Std. Error",
                                "t value", "Pr(>|t|)")]
  write.csv(coef_table, here(file_name), row.names = FALSE)
  cat("Saved:", file_name, "\n")
}

export_glm_output(baseline$model_fit, "baseline_glm_output.csv")
export_glm_output(stepwise$model_fit, "stepwise_glm_output.csv")
export_glm_output(boruta$model_fit,   "boruta_glm_output.csv")
export_glm_output(enet$model_fit,     "elastic_net_glm_output.csv")
