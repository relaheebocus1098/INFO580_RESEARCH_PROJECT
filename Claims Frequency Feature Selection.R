# =============================================================================
# CLAIMS FREQUENCY — FEATURE SELECTION COMPARISON
# Baseline NB GLM  |  Stepwise AIC  |  Boruta  |  Elastic Net
#
# Preprocessing mirrors the baseline modelling script:
#   - Exposure capped at 1 year (pmin)
#   - Client-stratified 80:20 train/test split (seed = 42)
#   - 5-fold grouped CV by client id (zero client leakage)
#   - Encoding recipe fitted on analysis partition only
#   - Premium features excluded; near-collinear history bands dropped
#   - Offset: log(exposure_years) in every NB GLM
#
# Selection criterion: NB deviance (coherent with MLE objective)
# Reporting metrics  : MAE, RMSE, Pseudo-R², NB Deviance — all as mean ± SE
# =============================================================================


# =============================================================================
# 0.  LIBRARIES
# =============================================================================

library(here)
library(dplyr)
library(rsample)       # initial_split, group_vfold_cv
library(recipes)       # recipe, step_*, prep, bake
library(MASS)          # glm.nb, stepAIC, negative.binomial family
library(glmnet)        # Elastic Net with NB family  (requires glmnet >= 4.0)
library(Boruta)        # Random Forest wrapper feature selection
library(gt)            # formatted HTML comparison table
library(tidyr)

# =============================================================================
# 1.  DATA LOADING & EXPOSURE CAPPING
# =============================================================================

basic_glm_data <- read.csv(
  here("Feature_Selection_Four_Wheels_GLM_R.csv")
)

# Exposure is capped at 1 year; all policies are practically annual.
basic_glm_data$exposure_years <- pmin(basic_glm_data$exposure_years, 1)


# =============================================================================
# 2.  VARIABLE TAXONOMY
# =============================================================================

# Nominal predictors — treatment-coded dummy expansion inside the recipe
one_hot_vars <- c(
  "area_label", "channel_label", "payment_label", "second_driver_label",
  "type_fuel", "seniority_category", "start_season", "start_day_of_week"
)

# Ordinal predictors — fixed level ordering applied before encoding
ord_vars <- c(
  "age_category", "age_when_licensed_band", "driving_exp_category",
  "vehicle_age_category", "power_band", "cylinder_band", "value_band",
  "weight_band", "length_imputed_band", "n_doors_band",
  # "claims_history_band"    # dropped: near-perfect collinearity
  # "r_claims_history_band"  # dropped: near-perfect collinearity
  "avg_claims_per_year_history_band", "cross_sell_ratio_band",
  "policy_retention_rate_band", "risk_score_band"
)

# Binary predictors — retained as 0/1 integers
binary_vars <- c(
  "growing_portfolio", "had_previous_lapse", "multi_policy_holder",
  "multi_product_holder", "start_weekend"
)


# =============================================================================
# 3.  CLIENT-STRATIFIED TRAIN / TEST SPLIT  (80:20)
# =============================================================================

# Split at client level — no client's policies appear in both partitions.
# Stratification is defined jointly by lapse history, claim history, and
# policy-count tier.

set.seed(42)

client_strat <- basic_glm_data %>%
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

# Collapse singleton strata so initial_split does not error
rare_strata  <- names(table(client_strat$strata)[table(client_strat$strata) == 1])
client_strat <- client_strat %>%
  mutate(strata = if_else(strata %in% rare_strata, "other", strata))

split_obj  <- initial_split(client_strat, prop = 0.8, strata = strata)
train_ids  <- training(split_obj)$id
test_ids   <- testing(split_obj)$id

train_data <- basic_glm_data %>% filter(id %in% train_ids)
test_data  <- basic_glm_data %>% filter(id %in% test_ids)

# Hard stop: zero client overlap between train and test
stopifnot(length(intersect(unique(train_data$id), unique(test_data$id))) == 0)
cat("Client leakage check (train/test): PASSED\n")
cat(sprintf("Train: %d policies | %d clients\n", nrow(train_data), n_distinct(train_data$id)))
cat(sprintf("Test : %d policies | %d clients\n", nrow(test_data),  n_distinct(test_data$id)))


# =============================================================================
# 4.  5-FOLD GROUPED CROSS-VALIDATION STRUCTURE
# =============================================================================

# Every policy belonging to a given client lands in the same fold.

set.seed(42)
cv_folds <- group_vfold_cv(train_data, group = id, v = 5)

leakage_free <- vapply(cv_folds$splits, function(sp) {
  length(intersect(unique(analysis(sp)$id), unique(assessment(sp)$id))) == 0L
}, logical(1))
cat(sprintf("CV client-leakage-free across all folds: %s\n", all(leakage_free)))


# =============================================================================
# 5.  ENCODING & EVALUATION HELPERS
# =============================================================================

# Fixed ordinal level orderings — consistent across folds and train/test
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

recode_ordinals <- function(df) {
  for (v in intersect(names(ordinal_levels), names(df)))
    df[[v]] <- factor(df[[v]], levels = ordinal_levels[[v]], ordered = FALSE)
  df
}

# Columns unconditionally excluded from every predictor set
always_exclude <- c(
  "id", "id2", "n_claims_year", "avg_claim_cost",
  "lapse_flag", "days_to_event", "exposure_years",
  "r_claims_history_band", "claims_history_band"
)

# Premium features withheld to prevent pricing-signal leakage
premium_exclude <- c("premium_band", "premium_vs_median_band")

get_predictor_cols <- function(data, target) {
  extra_drop <- if (target %in% c("n_claims_year", "avg_claim_cost"))
    premium_exclude else character(0)
  setdiff(names(data), unique(c(always_exclude, extra_drop, target)))
}

# Encoding pipeline fitted on the analysis partition only
encode_predictors <- function(train_df, test_df, target) {
  pred_cols <- get_predictor_cols(train_df, target)
  x_train   <- recode_ordinals(train_df[, pred_cols, drop = FALSE])
  x_test    <- recode_ordinals(test_df[,  pred_cols, drop = FALSE])

  rec <- recipe(~ ., data = x_train) %>%
    step_novel(all_nominal_predictors()) %>%
    step_unknown(all_nominal_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_zv(all_predictors())

  prep_rec <- prep(rec, training = x_train, retain = TRUE)
  list(
    train  = bake(prep_rec, new_data = NULL),
    test   = bake(prep_rec, new_data = x_test),
    recipe = prep_rec
  )
}

# Scalar metric functions
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae  <- function(obs, pred) mean(abs(obs  - pred),    na.rm = TRUE)

# Pseudo-R² = 1 - (residual deviance / null deviance)
pseudo_r2_glm <- function(fit) 1 - fit$deviance / fit$null.deviance

# Standard error of the mean across CV folds
cv_se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# NB log-likelihood on held-out predictions — used for out-of-sample deviance
# consistently across all methods, including Elastic Net.
nb_loglik <- function(y, mu, theta) {
  sum(
    lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
      y     * log(mu    / (mu + theta)) +
      theta * log(theta / (mu + theta))
  )
}

# Mean -2 × NB log-likelihood per observation. Lower is better.
nb_deviance_cv <- function(y, mu, theta) {
  -2 * nb_loglik(y, mu, theta) / length(y)
}

# Single-line fold progress printer
fold_done <- function(i, n_folds, method, mae_val, n_vars_val) {
  cat(sprintf("  [%d/%d] %-14s  MAE = %.4f  |  vars = %d\n",
              i, n_folds, method, mae_val, as.integer(round(n_vars_val))))
}

# Aggregate per-fold data frames into one mean ± SE summary row
build_cv_row <- function(cv_df, model_name) {
  if (is.list(cv_df) && !is.data.frame(cv_df)) cv_df <- bind_rows(cv_df)
  data.frame(
    Model       = model_name,
    CV_R2       = mean(cv_df$R2,       na.rm = TRUE),
    CV_R2_SE    = cv_se(cv_df$R2),
    CV_RMSE     = mean(cv_df$RMSE,     na.rm = TRUE),
    CV_RMSE_SE  = cv_se(cv_df$RMSE),
    CV_MAE      = mean(cv_df$MAE,      na.rm = TRUE),
    CV_MAE_SE   = cv_se(cv_df$MAE),
    CV_Dev      = mean(cv_df$Deviance, na.rm = TRUE),
    CV_Dev_SE   = cv_se(cv_df$Deviance),
    CV_Nvars    = mean(cv_df$n_vars,   na.rm = TRUE),
    CV_Nvars_SE = cv_se(cv_df$n_vars)
  )
}


# =============================================================================
# 6.  SHARED ENCODED MATRICES  (full train / test)
# =============================================================================

enc_all <- encode_predictors(train_data, test_data, "n_claims_year")
X_train <- enc_all$train
X_test  <- enc_all$test

train_mod <- cbind(
  y              = train_data$n_claims_year,
  exposure_years = train_data$exposure_years,
  X_train
)

test_mod <- cbind(
  y              = test_data$n_claims_year,
  exposure_years = test_data$exposure_years,
  X_test
)

oos_r2 <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}

# =============================================================================
# 7.  BASELINE: FULL NEGATIVE BINOMIAL GLM
# =============================================================================
# Full NB GLM on all available predictors; serves as the performance reference.
# theta_hat is passed as the fixed dispersion parameter to glmnet.

cat("\n====================================================\n")
cat("BASELINE: Full Negative Binomial GLM\n")
cat("====================================================\n")

nb_fit <- glm.nb(
  y ~ . - exposure_years + offset(log(exposure_years)),
  data = train_mod
)

theta_hat <- nb_fit$theta
cat(sprintf("Estimated dispersion theta: %.4f (SE: %.4f)\n",
            theta_hat, nb_fit$SE.theta))
cat(sprintf("Pearson dispersion ratio  : %.4f\n",
            sum(residuals(nb_fit, type = "pearson")^2) / df.residual(nb_fit)))

# --- CV: Baseline NB ---------------------------------------------------------
cat("  Running 5-fold CV...\n")
n_folds     <- length(cv_folds$splits)
cv_baseline <- lapply(seq_len(n_folds), function(i) {
  sp  <- cv_folds$splits[[i]]
  tr  <- analysis(sp);   va  <- assessment(sp)
  enc <- encode_predictors(tr, va, "n_claims_year")

  tr_mod <- cbind(y = tr$n_claims_year, exposure_years = tr$exposure_years, enc$train)
  va_mod <- cbind(y = va$n_claims_year, exposure_years = va$exposure_years, enc$test)

  m    <- glm.nb(y ~ . - exposure_years + offset(log(exposure_years)), data = tr_mod)
  pred <- predict(m, va_mod, type = "response")

  res <- data.frame(
    R2       = pseudo_r2_glm(m),
    RMSE     = rmse(va_mod$y, pred),
    MAE      = mae(va_mod$y,  pred),
    Deviance = nb_deviance_cv(va_mod$y, pred, m$theta),
    n_vars   = length(coef(m)) - 1L    # intercept excluded
  )
  fold_done(i, n_folds, "Baseline NB", res$MAE, res$n_vars)
  res
})

cv_baseline_df <- bind_rows(cv_baseline)

# Save fold-level results immediately — both formats
saveRDS(cv_baseline_df, here("cv_baseline_fold_results.rds"))
write.csv(cv_baseline_df, here("cv_baseline_fold_results.csv"), row.names = FALSE)
cat("  Saved: cv_baseline_fold_results  [.rds + .csv]\n")

cat("\nBaseline NB — CV results (mean ± SE):\n")
for (nm in c("R2", "RMSE", "MAE", "Deviance", "n_vars"))
  cat(sprintf("  %-10s  %.4f  ±  %.4f\n",
              nm,
              mean(cv_baseline_df[[nm]], na.rm = TRUE),
              cv_se(cv_baseline_df[[nm]])))

# --- Test set: Baseline ------------------------------------------------------
pred_nb_test  <- predict(nb_fit, test_mod, type = "response")

test_baseline <- data.frame(
  Model     = "Baseline NB",
  Test_MAE  = mae(test_mod$y,  pred_nb_test),
  Test_RMSE = rmse(test_mod$y, pred_nb_test),
  Test_Dev  = nb_deviance_cv(test_mod$y, pred_nb_test, theta_hat),
  Test_R2   = oos_r2(test_mod$y, pred_nb_test)
)

# Save the fitted baseline model object
saveRDS(nb_fit, here("model_baseline_nb.rds"))
cat("  Saved: model_baseline_nb.rds\n")


# =============================================================================
# 8.  FEATURE SELECTION METHOD 1: STEPWISE AIC
# =============================================================================
# stepAIC() performs bidirectional stepwise search over the NB GLM guided by
# AIC.  The full NB is fitted first on the analysis partition, then stepwise
# is applied within the same fold — the selected variable set is not informed
# by the validation partition in any fold.  The final production variable set
# is derived from the full-training-set run.

cat("\n====================================================\n")
cat("METHOD 1: Stepwise AIC\n")
cat("====================================================\n")

# --- CV: Stepwise AIC --------------------------------------------------------
cat("  Running 5-fold CV...\n")
cv_step <- lapply(seq_len(n_folds), function(i) {
  sp  <- cv_folds$splits[[i]]
  tr  <- analysis(sp);   va  <- assessment(sp)
  enc <- encode_predictors(tr, va, "n_claims_year")

  tr_mod <- cbind(y = tr$n_claims_year, exposure_years = tr$exposure_years, enc$train)
  va_mod <- cbind(y = va$n_claims_year, exposure_years = va$exposure_years, enc$test)

  # Fit full NB on analysis fold, then apply bidirectional stepwise
  full_m <- glm.nb(y ~ . - exposure_years + offset(log(exposure_years)), data = tr_mod)
  step_m <- stepAIC(full_m, direction = "both", trace = FALSE)

  pred <- predict(step_m, va_mod, type = "response")

  res <- data.frame(
    R2       = pseudo_r2_glm(step_m),
    RMSE     = rmse(va_mod$y, pred),
    MAE      = mae(va_mod$y,  pred),
    Deviance = nb_deviance_cv(va_mod$y, pred, step_m$theta),
    n_vars   = length(coef(step_m)) - 1L
  )
  fold_done(i, n_folds, "Stepwise AIC", res$MAE, res$n_vars)
  res
})

cv_step_df <- bind_rows(cv_step)

# Save fold-level results immediately — both formats
saveRDS(cv_step_df, here("cv_stepwise_fold_results.rds"))
write.csv(cv_step_df, here("cv_stepwise_fold_results.csv"), row.names = FALSE)
cat("  Saved: cv_stepwise_fold_results  [.rds + .csv]\n")

cat("\nStepwise AIC — CV results (mean ± SE):\n")
for (nm in c("R2", "RMSE", "MAE", "Deviance", "n_vars"))
  cat(sprintf("  %-10s  %.4f  ±  %.4f\n",
              nm,
              mean(cv_step_df[[nm]], na.rm = TRUE),
              cv_se(cv_step_df[[nm]])))

# --- Final model: Stepwise AIC on full training set --------------------------
nb_full_for_step <- glm.nb(
  y ~ . - exposure_years + offset(log(exposure_years)),
  data = train_mod
)
step_fit <- stepAIC(nb_full_for_step, direction = "both", trace = FALSE)

step_selected <- names(coef(step_fit))
step_selected <- step_selected[step_selected != "(Intercept)"]
cat(sprintf("\nStepwise AIC retained %d variables on full training set.\n",
            length(step_selected)))

# --- Test set: Stepwise AIC --------------------------------------------------
pred_step_test <- predict(step_fit, test_mod, type = "response")

test_step <- data.frame(
  Model     = "Stepwise AIC",
  Test_MAE  = mae(test_mod$y,  pred_step_test),
  Test_RMSE = rmse(test_mod$y, pred_step_test),
  Test_Dev  = nb_deviance_cv(test_mod$y, pred_step_test, step_fit$theta)
)

# Save the fitted stepwise model and its selected variable list
saveRDS(step_fit,      here("model_stepwise_nb.rds"))
saveRDS(step_selected, here("model_stepwise_selected_vars.rds"))
cat("  Saved: model_stepwise_nb.rds  |  model_stepwise_selected_vars.rds\n")


# =============================================================================
# 9.  FEATURE SELECTION METHOD 2: BORUTA
# =============================================================================
# Boruta is a Random Forest wrapper that iteratively compares each feature's
# importance to randomly shuffled shadow copies, confirming or rejecting each.
# TentativeRoughFix() resolves remaining tentative features using median
# importance.  The NB GLM is then refit on confirmed features with the offset.

cat("\n====================================================\n")
cat("METHOD 2: Boruta\n")
cat("====================================================\n")

# --- CV: Boruta --------------------------------------------------------------
cat("  Running 5-fold CV...\n")
cv_boruta <- lapply(seq_along(cv_folds$splits), function(i) {
  sp  <- cv_folds$splits[[i]]
  tr  <- analysis(sp);   va  <- assessment(sp)
  enc <- encode_predictors(tr, va, "n_claims_year")

  tr_mod <- cbind(y = tr$n_claims_year, exposure_years = tr$exposure_years, enc$train)
  va_mod <- cbind(y = va$n_claims_year, exposure_years = va$exposure_years, enc$test)

  set.seed(42 + i)
  bor <- Boruta(
    x       = enc$train,
    y       = tr$n_claims_year,
    maxRuns = 150,
    doTrace = 2,
    verbose = FALSE
  )
  
  bor       <- TentativeRoughFix(bor)
  confirmed <- getSelectedAttributes(bor, withTentative = FALSE)
  if (length(confirmed) == 0) confirmed <- names(enc$train)
  
  cat(sprintf("   -> Fold %d/%d completed — Refitting NB on %d confirmed vars...\n",
              i, n_folds, length(confirmed)))
  
  tr_sel <- tr_mod[, c("y", "exposure_years", confirmed), drop = FALSE]
  va_sel <- va_mod[, c("y", "exposure_years", confirmed), drop = FALSE]
  
  m    <- glm.nb(y ~ . - exposure_years + offset(log(exposure_years)), data = tr_sel)
  pred <- predict(m, va_sel, type = "response")
  
  res <- data.frame(
    R2       = pseudo_r2_glm(m),
    RMSE     = rmse(va_mod$y, pred),
    MAE      = mae(va_mod$y,  pred),
    Deviance = nb_deviance_cv(va_mod$y, pred, m$theta),
    n_vars   = length(confirmed)
  )
  fold_done(i, n_folds, "Boruta", res$MAE, res$n_vars)
  res
})

cat("\n  Running final full Boruta on training set...\n")
set.seed(42)
bor_full <- Boruta(
  x       = X_train,
  y       = train_data$n_claims_year,
  maxRuns = 150,
  doTrace = 2,
  verbose = FALSE
)

bor_full        <- TentativeRoughFix(bor_full)
boruta_selected <- getSelectedAttributes(bor_full, withTentative = FALSE)
cat(sprintf("\nBoruta confirmed %d variables on full training set.\n",
            length(boruta_selected)))
train_bor <- train_mod[, c("y", "exposure_years", boruta_selected), drop = FALSE]
test_bor  = test_mod[,  c("y", "exposure_years", boruta_selected), drop = FALSE]

boruta_fit <- glm.nb(
  y ~ . - exposure_years + offset(log(exposure_years)),
  data = train_bor
)

# --- Test set: Boruta --------------------------------------------------------
pred_boruta_test <- predict(boruta_fit, test_bor, type = "response")

test_boruta_res <- data.frame(
  Model     = "Boruta",
  Test_MAE  = mae(test_mod$y,  pred_boruta_test),
  Test_RMSE = rmse(test_mod$y, pred_boruta_test),
  Test_Dev  = nb_deviance_cv(test_mod$y, pred_boruta_test, boruta_fit$theta)
)

# Save the fitted Boruta NB model, the Boruta object, and the selected variable list
saveRDS(boruta_fit,      here("model_boruta_nb.rds"))
saveRDS(bor_full,        here("model_boruta_object.rds"))
saveRDS(boruta_selected, here("model_boruta_selected_vars.rds"))
cat("  Saved: model_boruta_nb.rds  |  model_boruta_object.rds  |  model_boruta_selected_vars.rds\n")

cv_boruta_df <- bind_rows(cv_boruta)
saveRDS(cv_boruta_df, here("cv_boruta_fold_results.rds"))
write.csv(cv_boruta_df, here("cv_boruta_fold_results.csv"), row.names = FALSE)
cat("  Saved: cv_boruta_fold_results  [.rds + .csv]\n")


# =============================================================================
# 10. FEATURE SELECTION METHOD 3: ELASTIC NET
# =============================================================================
# glmnet fits a penalised NB GLM via coordinate descent.
# alpha controls the Lasso/Ridge mix; lambda controls regularisation strength
# and is selected by internal 5-fold CV using NB deviance. lambda.1se is used
# for parsimony over lambda.min.  theta_hat from the baseline is passed as
# the fixed dispersion parameter.  The outer loop iterates over 5 grouped CV
# folds; best alpha is the value minimising mean outer-fold CV deviance.

cat("\n====================================================\n")
cat("METHOD 3: Elastic Net\n")
cat("====================================================\n")

alpha_grid <- c(0.1, 0.25, 0.5, 0.75, 1.0)

# --- CV: Elastic Net (outer grouped folds × alpha grid) ----------------------
# One progress line per fold; alpha sub-steps are silent to avoid crowding.
cat(sprintf("  Running 5-fold CV over alpha grid {%s}...\n",
            paste(alpha_grid, collapse = ", ")))

cv_enet_all <- lapply(seq_len(n_folds), function(i) {
  sp  <- cv_folds$splits[[i]]
  tr  <- analysis(sp);   va  <- assessment(sp)
  enc <- encode_predictors(tr, va, "n_claims_year")

  X_tr   <- as.matrix(enc$train)
  X_va   <- as.matrix(enc$test)
  y_tr   <- tr$n_claims_year
  y_va   <- va$n_claims_year
  off_tr <- log(tr$exposure_years)
  off_va <- log(va$exposure_years)

  # For each alpha: select lambda via internal cv.glmnet on analysis partition
  fold_results <- lapply(alpha_grid, function(a) {
    set.seed(42)
    cv_fit <- cv.glmnet(
      x            = X_tr,
      y            = y_tr,
      family       = MASS::negative.binomial(theta = theta_hat),
      alpha        = a,
      offset       = off_tr,
      type.measure = "deviance",
      nfolds       = 5
    )

    pred <- as.vector(
      predict(cv_fit, newx = X_va, newoffset = off_va,
              s = "lambda.1se", type = "response")
    )

    coefs     <- coef(cv_fit, s = "lambda.1se")
    n_nonzero <- sum(coefs[-1] != 0)

    data.frame(
      alpha    = a,
      R2       = 1 - sum((y_va - pred)^2) / sum((y_va - mean(y_va))^2),
      RMSE     = rmse(y_va, pred),
      MAE      = mae(y_va,  pred),
      Deviance = nb_deviance_cv(y_va, pred, theta_hat),
      n_vars   = n_nonzero
    )
  })

  res <- bind_rows(fold_results)
  best_row <- res[which.min(res$Deviance), ]
  cat(sprintf("  [%d/%d] Elastic Net       MAE = %.4f  |  best alpha = %.2f  |  vars = %d\n",
              i, n_folds, best_row$MAE, best_row$alpha, as.integer(best_row$n_vars)))
  res
})

cv_enet_df <- bind_rows(cv_enet_all)

# Save fold-level results immediately — both formats
saveRDS(cv_enet_df, here("cv_elasticnet_fold_results.rds"))
write.csv(cv_enet_df, here("cv_elasticnet_fold_results.csv"), row.names = FALSE)
cat("  Saved: cv_elasticnet_fold_results  [.rds + .csv]\n")

# Summarise CV metrics by alpha value across the 5 outer folds
enet_alpha_summary <- cv_enet_df %>%
  group_by(alpha) %>%
  summarise(
    MAE_mean  = mean(MAE),       MAE_se   = cv_se(MAE),
    RMSE_mean = mean(RMSE),      RMSE_se  = cv_se(RMSE),
    Dev_mean  = mean(Deviance),  Dev_se   = cv_se(Deviance),
    R2_mean   = mean(R2),        R2_se    = cv_se(R2),
    nvar_mean = mean(n_vars),    nvar_se  = cv_se(n_vars),
    .groups = "drop"
  )

cat("\nElastic Net — CV summary by alpha (mean ± SE across 5 outer folds):\n")
print(as.data.frame(enet_alpha_summary), digits = 4)

# Select the alpha minimising mean outer-fold CV deviance
best_alpha <- enet_alpha_summary$alpha[which.min(enet_alpha_summary$Dev_mean)]
cat(sprintf("\nSelected alpha: %.2f  (minimum CV deviance)\n", best_alpha))

# Extract per-fold rows for the best alpha to feed into the summary table
cv_enet_best <- cv_enet_df %>% filter(alpha == best_alpha)

# --- Final model: Elastic Net on full training set ---------------------------
X_tr_full   <- as.matrix(X_train)
X_te_full   <- as.matrix(X_test)
off_tr_full <- log(train_data$exposure_years)
off_te_full <- log(test_data$exposure_years)

set.seed(42)
cv_enet_final <- cv.glmnet(
  x            = X_tr_full,
  y            = train_data$n_claims_year,
  family       = MASS::negative.binomial(theta = theta_hat),
  alpha        = best_alpha,
  offset       = off_tr_full,
  type.measure = "deviance",
  nfolds       = 5
)

enet_coefs         <- coef(cv_enet_final, s = "lambda.1se")
enet_selected      <- rownames(enet_coefs)[enet_coefs[, 1] != 0]
enet_selected      <- enet_selected[enet_selected != "(Intercept)"]
cat(sprintf("Elastic Net retained %d variables at lambda.1se (alpha = %.2f).\n",
            length(enet_selected), best_alpha))

# --- Test set: Elastic Net ---------------------------------------------------
pred_enet_test <- as.vector(
  predict(cv_enet_final, newx = X_te_full, newoffset = off_te_full,
          s = "lambda.1se", type = "response")
)

test_enet <- data.frame(
  Model     = "Elastic Net",
  Test_MAE  = mae(test_mod$y,  pred_enet_test),
  Test_RMSE = rmse(test_mod$y, pred_enet_test),
  Test_Dev  = nb_deviance_cv(test_mod$y, pred_enet_test, theta_hat)
)

# Save the fitted Elastic Net cv.glmnet object and the selected variable list
saveRDS(cv_enet_final, here("model_elasticnet_cvglmnet.rds"))
saveRDS(enet_selected, here("model_elasticnet_selected_vars.rds"))
cat(sprintf("  Saved: model_elasticnet_cvglmnet.rds  |  model_elasticnet_selected_vars.rds\n"))


# =============================================================================
# RELOAD SAVED OBJECTS AND RECONSTRUCT TEST ROWS
# =============================================================================

cv_baseline_df <- readRDS(here("cv_baseline_fold_results.rds"))
cv_step_df     <- readRDS(here("cv_stepwise_fold_results.rds"))
cv_boruta_df   <- readRDS(here("cv_boruta_fold_results.rds"))
cv_enet_best   <- readRDS(here("cv_elasticnet_fold_results.rds"))

# Reconstruct test rows from saved models
# 1. Baseline
baseline_fit  <- readRDS(here("model_baseline_nb.rds"))
pred_baseline <- predict(baseline_fit, test_mod, type = "response")
test_baseline <- data.frame(
  Model     = "Baseline NB",
  Test_MAE  = mae(test_mod$y,  pred_baseline),
  Test_RMSE = rmse(test_mod$y, pred_baseline),
  Test_Dev  = nb_deviance_cv(test_mod$y, pred_baseline, baseline_fit$theta),
  Test_R2   = oos_r2(test_mod$y, pred_baseline)
)

# 2. Stepwise AIC
stepwise_fit  <- readRDS(here("model_stepwise_nb.rds"))
step_selected <- readRDS(here("model_stepwise_selected_vars.rds"))
test_step_df  <- test_mod[, c("y", "exposure_years", step_selected), drop = FALSE]
pred_step     <- predict(stepwise_fit, test_step_df, type = "response")
test_step     <- data.frame(
  Model     = "Stepwise AIC",
  Test_MAE  = mae(test_mod$y,  pred_step),
  Test_RMSE = rmse(test_mod$y, pred_step),
  Test_Dev  = nb_deviance_cv(test_mod$y, pred_step, stepwise_fit$theta),
  Test_R2   = oos_r2(test_mod$y, pred_step)
)

# 3. Boruta
boruta_fit      <- readRDS(here("model_boruta_nb.rds"))
boruta_selected <- readRDS(here("model_boruta_selected_vars.rds"))
test_bor        <- test_mod[, c("y", "exposure_years", boruta_selected), drop = FALSE]
pred_boruta     <- predict(boruta_fit, test_bor, type = "response")
test_boruta_res <- data.frame(
  Model     = "Boruta",
  Test_MAE  = mae(test_mod$y,  pred_boruta),
  Test_RMSE = rmse(test_mod$y, pred_boruta),
  Test_Dev  = nb_deviance_cv(test_mod$y, pred_boruta, boruta_fit$theta),
  Test_R2   = oos_r2(test_mod$y, pred_boruta)
)

# 4. Elastic Net
enet_fit      <- readRDS(here("model_elasticnet_cvglmnet.rds"))
enet_selected <- readRDS(here("model_elasticnet_selected_vars.rds"))

enet_vars  <- rownames(coef(enet_fit))[-1]
X_test_mat <- data.matrix(X_test[, enet_vars, drop = FALSE])

if (!is.null(enet_fit$glmnet.fit$offset) && enet_fit$glmnet.fit$offset) {
  pred_enet <- predict(enet_fit, newx = X_test_mat, s = "lambda.1se",
                       newoffset = log(test_mod$exposure_years), type = "response")[, 1]
} else {
  pred_enet <- predict(enet_fit, newx = X_test_mat, s = "lambda.1se", type = "response")[, 1]
}

t_theta <- theta_hat
test_enet <- data.frame(
  Model     = "Elastic Net",
  Test_MAE  = mae(test_mod$y,  pred_enet),
  Test_RMSE = rmse(test_mod$y, pred_enet),
  Test_Dev  = nb_deviance_cv(test_mod$y, pred_enet, t_theta),
  Test_R2   = oos_r2(test_mod$y, pred_enet)
)




# =============================================================================
# 11. ASSEMBLE FULL RESULTS TABLE
# =============================================================================

cv_summary <- bind_rows(
  build_cv_row(cv_baseline_df, "Baseline NB"),
  build_cv_row(cv_step_df,     "Stepwise AIC"),
  build_cv_row(cv_boruta_df,   "Boruta"),
  build_cv_row(cv_enet_best,   "Elastic Net")
)

test_summary <- bind_rows(test_baseline, test_step, test_boruta_res, test_enet)

results_full <- left_join(cv_summary, test_summary, by = "Model")

cat("\n====================================================\n")
cat("FULL RESULTS SUMMARY\n")
cat("====================================================\n")
print(results_full, digits = 4)


# =============================================================================
# 12. EXPORT TO CSV
# =============================================================================
write.csv(results_full, here("feature_selection_cv_results.csv"), row.names = FALSE)
cat("\nExported: feature_selection_cv_results.csv\n")

if (exists("enet_alpha_summary")) {
  write.csv(enet_alpha_summary, here("elastic_net_alpha_cv_summary.csv"), row.names = FALSE)
  cat("Exported: elastic_net_alpha_cv_summary.csv\n")
} else {
  cat("Note: enet_alpha_summary row details not in current environment; skipped summary csv.\n")
}

# =============================================================================
# 13. FINAL COMPARISON TABLE  (gt — saves as HTML)
# =============================================================================

fmt_pm <- function(m, s, d = 4) sprintf(paste0("%.", d, "f \u00b1 %.", d, "f"), m, s)

best_cv_r2_model   <- results_full$Model[which.max(results_full$CV_R2)]
best_cv_rmse_model <- results_full$Model[which.min(results_full$CV_RMSE)]
best_cv_mae_model  <- results_full$Model[which.min(results_full$CV_MAE)]
best_cv_dev_model  <- results_full$Model[which.min(results_full$CV_Dev)]

best_test_mae_model  <- results_full$Model[which.min(results_full$Test_MAE)]
best_test_rmse_model <- results_full$Model[which.min(results_full$Test_RMSE)]
best_test_dev_model  <- results_full$Model[which.min(results_full$Test_Dev)]
best_test_r2_model   <- results_full$Model[which.max(results_full$Test_R2)]

display_tbl <- results_full %>%
  mutate(
    CV_R2_fmt    = fmt_pm(CV_R2,    CV_R2_SE),
    CV_RMSE_fmt  = fmt_pm(CV_RMSE,  CV_RMSE_SE),
    CV_MAE_fmt   = fmt_pm(CV_MAE,   CV_MAE_SE),
    CV_Dev_fmt   = fmt_pm(CV_Dev,   CV_Dev_SE),
    CV_Nvars_fmt = fmt_pm(CV_Nvars, CV_Nvars_SE, d = 1),
    Test_MAE_fmt  = round(Test_MAE,  4),
    Test_RMSE_fmt = round(Test_RMSE, 4),
    Test_Dev_fmt  = round(Test_Dev,  4),
    Test_R2_fmt   = round(Test_R2,   4) # The newly added Test R2
  ) %>%
  dplyr::select(
    Model,
    CV_R2_fmt, CV_RMSE_fmt, CV_MAE_fmt, CV_Dev_fmt, CV_Nvars_fmt,
    Test_MAE_fmt, Test_RMSE_fmt, Test_Dev_fmt, Test_R2_fmt
  )

# Build the beautiful gt table
# Build the beautiful gt table
gt_table <- display_tbl %>%
  gt(rowname_col = "Model") %>%
  
  # Assign the pretty Unicode/Markdown labels safely here
  cols_label(
    CV_R2_fmt    = md("CV R&sup2;"),
    CV_RMSE_fmt  = md("CV RMSE"),
    CV_MAE_fmt   = md("CV MAE"),
    CV_Dev_fmt   = md("CV Deviance"),
    CV_Nvars_fmt = md("CV N Vars"),
    Test_MAE_fmt  = md("Test MAE"),
    Test_RMSE_fmt = md("Test RMSE"),
    Test_Dev_fmt  = md("Test Deviance"),
    Test_R2_fmt   = md("Test R&sup2;")
  ) %>%
  
  # Title & subtitle
  tab_header(
    title    = md("**Claims Frequency — Feature Selection Comparison**"),
    subtitle = md(
      "Negative Binomial GLM &nbsp;|&nbsp; 5-Fold Grouped CV &nbsp;|&nbsp; 
       Client-stratified 80:20 split &nbsp;|&nbsp; seed = 42"
    )
  ) %>%
  
  # Column spanners
  tab_spanner(
    label   = md("**Cross-Validation (5-fold, mean \u00b1 SE)**"),
    columns = c(CV_R2_fmt, CV_RMSE_fmt, CV_MAE_fmt, CV_Dev_fmt, CV_Nvars_fmt)
  ) %>%
  tab_spanner(
    label   = md("**Hold-Out Test Set**"),
    columns = c(Test_MAE_fmt, Test_RMSE_fmt, Test_Dev_fmt, Test_R2_fmt)
  ) %>%
  
  # Bold row labels
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_stub()
  ) %>%
  
  # Shade Baseline NB as reference row background
  tab_style(
    style     = cell_fill(color = "#D6EAF8"),
    locations = cells_body(rows = Model == "Baseline NB")
  ) %>%
  
  # ---------------------------------------------------------------------------
  # Highlight best metric in each column
  tab_style(
  style     = cell_fill(color = "#D5F5E3"),
  locations = cells_body(columns = CV_R2_fmt, rows = Model == best_cv_r2_model)
) %>%
  tab_style(
    style     = cell_fill(color = "#D5F5E3"),
    locations = cells_body(columns = CV_RMSE_fmt, rows = Model == best_cv_rmse_model)
  ) %>%
  tab_style(
    style     = cell_fill(color = "#D5F5E3"),
    locations = cells_body(columns = CV_MAE_fmt, rows = Model == best_cv_mae_model)
  ) %>%
  tab_style(
    style     = cell_fill(color = "#D5F5E3"),
    locations = cells_body(columns = CV_Dev_fmt, rows = Model == best_cv_dev_model)
  ) %>%
  
  # Test Highlights
  tab_style(
    style     = cell_fill(color = "#D5F5E3"),
    locations = cells_body(columns = Test_MAE_fmt, rows = Model == best_test_mae_model)
  ) %>%
  tab_style(
    style     = cell_fill(color = "#D5F5E3"),
    locations = cells_body(columns = Test_RMSE_fmt, rows = Model == best_test_rmse_model)
  ) %>%
  tab_style(
    style     = cell_fill(color = "#D5F5E3"),
    locations = cells_body(columns = Test_Dev_fmt, rows = Model == best_test_dev_model)
  ) %>%
  tab_style(
    style     = cell_fill(color = "#D5F5E3"),
    locations = cells_body(columns = Test_R2_fmt, rows = Model == best_test_r2_model)
  ) %>%


# Centre all value columns
cols_align(align = "center", columns = everything()) %>%
  
  # Footnotes
  tab_footnote(
    footnote  = md(
      "CV Deviance = mean &minus;2 &times; NB log-likelihood per observation 
       on the validation fold. Lower is better."
    ),
    locations = cells_column_labels(columns = CV_Dev_fmt)
  ) %>%
  tab_footnote(
    footnote  = md(
      paste0(
        "Elastic Net: alpha selected from grid {",
        paste(alpha_grid, collapse = ", "),
        "} by minimum CV deviance; lambda selected via lambda.1se."
      )
    ),
    locations = cells_stub(rows = "Elastic Net")
  ) %>%
  tab_footnote(
    footnote  = md(
      "Green cells indicate the mathematically ideal performance metric in that column."
    ),
    locations = cells_column_labels(columns = CV_MAE_fmt)
  ) %>%
  
  # Styling
  tab_options(
    table.font.size                 = px(13),
    heading.title.font.size         = px(15),
    heading.subtitle.font.size      = px(12),
    column_labels.font.weight       = "bold",
    stub.font.weight                = "bold",
    table.border.top.color          = "#2C3E50",
    table.border.top.width          = px(2),
    table.border.bottom.color       = "#2C3E50",
    table.border.bottom.width       = px(2),
    heading.border.bottom.color     = "#2C3E50",
    row.striping.include_table_body = TRUE,
    row.striping.background_color   = "#F8F9FA"
  )

gtsave(gt_table, here("feature_selection_comparison_table.html"))
cat("Saved: feature_selection_comparison_table.html\n")

# Display in RStudio Viewer pane
gt_table

# =============================================================================
# 14. FINAL STEPWISE MODEL SUMMARY
# =============================================================================
# Refit the selected Stepwise AIC model on the full training set and export
# the coefficient table. This is the production frequency model.

step_selected  <- readRDS(here("model_stepwise_selected_vars.rds"))

formula_str   <- paste("y ~", paste(step_selected, collapse = " + "),
                       "+ offset(log(exposure_years))")
final_formula <- as.formula(formula_str)

final_nb_model <- glm.nb(final_formula, data = train_mod)
summary(final_nb_model)

coef_table          <- as.data.frame(summary(final_nb_model)$coefficients)
coef_table$Variable <- rownames(coef_table)
coef_table          <- coef_table[, c("Variable", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
colnames(coef_table) <- c("Variable", "Estimate", "Std. Error", "z value", "Pr(>|z|)")

write.csv(coef_table, here("stepwise_model_coefficients.csv"), row.names = FALSE)
cat("Exported: stepwise_model_coefficients.csv\n")

