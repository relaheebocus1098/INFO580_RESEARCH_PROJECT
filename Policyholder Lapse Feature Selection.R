# ==============================================================================
# LAPSE MODELLING — FEATURE SELECTION
# Methods : Stepwise AIC | Boruta | Elastic Net
# Models  : Stage 1 (lapse_flag) | Cox PH (days_to_event) | Stage 2 (early_lapse)
# Training metric (primary)  : AUC  (binary models) | C-Index (Cox)
# Training metric (secondary): Brier Score (binary) | IBS     (Cox)
# All CV metrics reported as Mean +/- SE across the same 5 client-stratified folds
# used throughout the baseline analysis.
#
# CV ARCHITECTURE (fold-blind selection — consistent with frequency/severity files)
# ---------------------------------------------------------------------------------
# Selection algorithm runs INSIDE each CV fold on the fold's own encoded data.
# Validation rows are never seen during that fold's feature selection.
#   Stepwise AIC : full GLM/Cox fit → stepAIC inside each fold
#   Boruta       : Boruta() inside each fold; doTrace=2, verbose=FALSE
#   Elastic Net  : outer fold × ALPHA_GRID for alpha selection;
#                  then cv.glmnet with best_alpha inside each fold for variable sets
# A separate full-training run of each algorithm produces the final feature
# set used for test-set evaluation and the saved production model.
# ==============================================================================

# ==============================================================================
# SECTION 0: LIBRARIES
# ==============================================================================

library(here)
library(dplyr)
library(rsample)
library(recipes)
library(MASS)        # stepAIC
library(survival)
library(pROC)
library(timeROC)
library(glmnet)      # Elastic Net
library(Boruta)      # Boruta feature selection
library(knitr)
library(kableExtra)  # HTML comparison tables

# ==============================================================================
# SECTION 1: UTILITY HELPERS
# ==============================================================================

progress <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush.console()
}

save_results <- function(obj, method, model_name, type = "metrics") {
  tag  <- sprintf("fs_%s_%s_%s", method, model_name, type)
  rds  <- here(paste0(tag, ".rds"))
  csv  <- here(paste0(tag, ".csv"))
  saveRDS(obj, rds)
  if (is.data.frame(obj)) write.csv(obj, csv, row.names = FALSE)
  progress(sprintf("  Saved: %s.rds / .csv", tag))
}

# ==============================================================================
# SECTION 2: METRIC HELPERS
# ==============================================================================

LAPSE_THRESHOLD <- 0.5
LANDMARK_T1     <- 300L
LANDMARK_T2     <- 365L
BRIER_GRID      <- c(seq(30L, 330L, by = 30L), 365L)

# suppressMessages() double-guards against pROC version differences;
# quiet=TRUE is the pROC >= 1.15 mechanism; suppressMessages covers older versions.
# Together they prevent "Setting levels: control = 0, case = 1" from printing.
auc_score <- function(obs, pred) {
  if (length(unique(obs)) < 2) return(NA_real_)
  suppressMessages(as.numeric(pROC::auc(obs, pred, quiet = TRUE)))
}

cindex_score <- function(time, event, lp)
  survival::concordance(Surv(time, event) ~ I(-lp))$concordance

rmse_prob      <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae_prob       <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
brier_score    <- function(obs, pred) mean((obs - pred)^2, na.rm = TRUE)
binom_deviance <- function(obs, pred) {
  pred <- pmax(pmin(pred, 1 - 1e-10), 1e-10)
  -2 * sum(obs * log(pred) + (1 - obs) * log(1 - pred), na.rm = TRUE)
}
corr_r2   <- function(obs, pred) cor(obs, pred, use = "complete.obs")^2
gini_coef <- function(auc) 2 * auc - 1

confusion_metrics <- function(obs, pred, threshold = LAPSE_THRESHOLD) {
  pc <- as.integer(pred >= threshold)
  TP <- sum(pc == 1 & obs == 1, na.rm = TRUE)
  TN <- sum(pc == 0 & obs == 0, na.rm = TRUE)
  FP <- sum(pc == 1 & obs == 0, na.rm = TRUE)
  FN <- sum(pc == 0 & obs == 1, na.rm = TRUE)
  N  <- TP + TN + FP + FN
  Pr <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  Re <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  Sp <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  F1 <- if (!is.na(Pr) && !is.na(Re) && (Pr + Re) > 0)
    2 * Pr * Re / (Pr + Re) else NA_real_
  data.frame(TP=TP,TN=TN,FP=FP,FN=FN,
             Accuracy=(TP+TN)/N,Precision=Pr,Recall=Re,Specificity=Sp,F1=F1)
}

binary_lapse_metrics <- function(obs, pred) {
  auc <- auc_score(obs, pred); cm <- confusion_metrics(obs, pred)
  data.frame(R2=corr_r2(obs,pred),RMSE=rmse_prob(obs,pred),
             MAE=mae_prob(obs,pred),Brier=brier_score(obs,pred),
             Binom_Dev=binom_deviance(obs,pred),
             AUC=auc,Gini=gini_coef(auc),
             Accuracy=cm$Accuracy,Precision=cm$Precision,
             Recall=cm$Recall,Specificity=cm$Specificity,F1=cm$F1)
}

cv_summary <- function(tbl) {
  means <- colMeans(tbl, na.rm = TRUE)
  sds   <- apply(tbl, 2, sd, na.rm = TRUE)
  n     <- colSums(!is.na(tbl))
  rbind(Mean = means, SE = sds / sqrt(n))
}

get_surv_from_cox <- function(cox_fit, lp_vec) {
  bh <- basehaz(cox_fit, centered = TRUE)
  H0 <- function(t) {
    idx <- which(bh$time <= t)
    if (!length(idx)) 0 else bh$hazard[max(idx)]
  }
  data.frame(S_T1 = exp(-H0(LANDMARK_T1) * exp(lp_vec)),
             S_T2 = exp(-H0(LANDMARK_T2) * exp(lp_vec)))
}

ipcw_brier_score <- function(train_time, train_event, val_time, val_event,
                             surv_pred, t_star) {
  km_cens <- survfit(Surv(train_time, 1L - train_event) ~ 1)
  G_fn    <- stepfun(km_cens$time, c(1, km_cens$surv))
  G_tstar <- G_fn(t_star)
  if (is.na(G_tstar) || G_tstar <= 0) return(NA_real_)
  n <- length(val_time); bsv <- rep(NA_real_, n)
  early <- val_event == 1L & val_time <= t_star
  if (any(early)) {
    idx <- which(early); G_ti <- G_fn(val_time[idx])
    ok  <- !is.na(G_ti) & G_ti > 0
    bsv[idx[ok]] <- (0 - surv_pred[idx[ok]])^2 / G_ti[ok]
  }
  past <- val_time > t_star
  if (any(past)) bsv[which(past)] <- (1 - surv_pred[which(past)])^2 / G_tstar
  mean(bsv, na.rm = TRUE)
}

integrated_brier_score <- function(train_time, train_event, val_time, val_event,
                                   cox_fit, lp_val, grid = BRIER_GRID) {
  bh <- basehaz(cox_fit, centered = TRUE)
  bs <- sapply(grid, function(t) {
    idx <- which(bh$time <= t)
    H0  <- if (!length(idx)) 0 else bh$hazard[max(idx)]
    ipcw_brier_score(train_time, train_event, val_time, val_event,
                     exp(-H0 * exp(lp_val)), t)
  })
  valid <- !is.na(bs)
  if (sum(valid) < 2) return(NA_real_)
  g <- grid[valid]; bv <- bs[valid]
  sum(diff(g) * (bv[-length(bv)] + bv[-1]) / 2) / (max(g) - min(g))
}

td_auc_at <- function(time_vec, event_vec, lp_vec, t_star) {
  tryCatch({
    res <- timeROC::timeROC(T=time_vec, delta=event_vec, marker=lp_vec,
                            cause=1L, times=t_star, iid=FALSE)
    av  <- as.numeric(res$AUC); av <- av[!is.na(av) & av > 0]
    if (!length(av)) NA_real_ else av[length(av)]
  }, error = function(e) NA_real_)
}

# ==============================================================================
# SECTION 3: ENCODING SCHEMA
# ==============================================================================

ordinal_levels <- list(
  age_category=c("18-25","26-35","36-45","46-55","56-65","65+"),
  age_when_licensed_band=c("<=25",">25"),
  driving_exp_category=c("<=5yrs","6-10yrs","11-20yrs","20+yrs"),
  vehicle_age_category=c("0-3yrs","4-7yrs","8-12yrs","12+yrs"),
  power_band=c("<100",">100"),
  cylinder_band=c("0-1400","1400-1600","1600-2000","2000+"),
  value_band=c("0-17500","17500-25000","25000-32500","32500+"),
  weight_band=c("0-1100","1100-1500","1500+"),
  length_imputed_band=c("<4","4-4.5",">4.5"),
  n_doors_band=c("<=4",">4"),
  claims_history_band=c("zero","one","two_to_three","four_to_five","six_plus"),
  r_claims_history_band=c("0","0-0.35","0.35-0.7","0.7-1.5","1.5+"),
  avg_claims_per_year_history_band=c("0","0-1",">1"),
  cross_sell_ratio_band=c("<0.75",">=0.75"),
  policy_retention_rate_band=c("<1","1"),
  risk_score_band=c("0","0-0.01","0.01-0.15","0.15-0.2","0.2+"),
  premium_band=c("0-250","250-300","300-350","350-400","400+"),
  premium_vs_median_band=c("-15%","-15% to 0%","0% to 15%","15% to 30%","30%+")
)

recode_ordinals <- function(df) {
  for (v in intersect(names(ordinal_levels), names(df)))
    df[[v]] <- factor(df[[v]], levels = ordinal_levels[[v]], ordered = FALSE)
  df
}

always_exclude <- c("id","id2","n_claims_year","avg_claim_cost_capped",
                    "lapse_flag","days_to_event","exposure_years",
                    "r_claims_history_band","claims_history_band")

get_predictor_cols <- function(data, target)
  setdiff(names(data), unique(c(always_exclude, target)))

encode_predictors <- function(train_df, test_df, target) {
  pred_cols <- get_predictor_cols(train_df, target)
  x_tr <- recode_ordinals(train_df[, pred_cols, drop = FALSE])
  x_te <- recode_ordinals(test_df[,  pred_cols, drop = FALSE])
  rec  <- recipe(~ ., data = x_tr) %>%
    step_novel(all_nominal_predictors()) %>%
    step_unknown(all_nominal_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_zv(all_predictors())
  prp <- prep(rec, training = x_tr, retain = TRUE)
  list(train = bake(prp, new_data = NULL),
       test  = bake(prp, new_data = x_te),
       recipe = prp)
}

strip_suffix <- function(dummy_names, known_vars = NULL) {
  if (is.null(known_vars)) return(sub("_[^_]+$", "", dummy_names))
  sapply(dummy_names, function(nm) {
    cands <- known_vars[startsWith(nm, paste0(known_vars, "_"))]
    if (length(cands) > 0) cands[which.max(nchar(cands))] else nm
  }, USE.NAMES = FALSE)
}

# ==============================================================================
# SECTION 4: DATA LOAD AND PREPROCESSING
# ==============================================================================

progress("Loading data and running preprocessing...")

basic_glm_data_with_avg_claim_cost <- read.csv(
  here("Feature_Selection_Four_Wheels_GLM_R.csv")
)
basic_glm_data_with_avg_claim_cost$avg_claim_cost_capped <-
  pmin(basic_glm_data_with_avg_claim_cost$avg_claim_cost, 1500)
basic_glm_data <- subset(basic_glm_data_with_avg_claim_cost, select = -avg_claim_cost)

set.seed(42)
client_strat <- basic_glm_data %>%
  group_by(id) %>%
  summarise(ever_lapsed  = as.integer(any(lapse_flag == 1, na.rm = TRUE)),
            ever_claimed = as.integer(any(n_claims_year > 0, na.rm = TRUE)),
            n_policies   = n(), .groups = "drop") %>%
  mutate(policy_tier = case_when(n_policies==1~"1",n_policies==2~"2",TRUE~"3+"),
         strata = paste(ever_lapsed, ever_claimed, policy_tier, sep="_"))
rare_strata  <- names(table(client_strat$strata)[table(client_strat$strata)==1])
client_strat <- client_strat %>%
  mutate(strata = if_else(strata %in% rare_strata, "other", strata))

split_obj  <- initial_split(client_strat, prop = 0.8, strata = strata)
train_data <- basic_glm_data %>% filter(id %in% training(split_obj)$id)
test_data  <- basic_glm_data %>% filter(id %in% testing(split_obj)$id)

set.seed(42)
cv_folds <- group_vfold_cv(train_data, group = id, v = 5)

train_data_lapsed <- train_data %>% filter(lapse_flag == 1) %>%
  mutate(early_lapse = as.integer(days_to_event <= LANDMARK_T1))
test_data_lapsed  <- test_data  %>% filter(lapse_flag == 1) %>%
  mutate(early_lapse = as.integer(days_to_event <= LANDMARK_T1))

enc_s1 <- encode_predictors(train_data,        test_data,        "lapse_flag")
enc_cx <- encode_predictors(train_data,        test_data,        "days_to_event")
enc_s2 <- encode_predictors(train_data_lapsed, test_data_lapsed, "early_lapse")

train_mod_lap <- cbind(y=train_data$lapse_flag,        enc_s1$train)
test_mod_lap  <- cbind(y=test_data$lapse_flag,         enc_s1$test)
train_mod_srv <- cbind(time=train_data$days_to_event,
                       event=train_data$lapse_flag,   enc_cx$train)
test_mod_srv  <- cbind(time=test_data$days_to_event,
                       event=test_data$lapse_flag,    enc_cx$test)
train_s2      <- cbind(y=train_data_lapsed$early_lapse, enc_s2$train)
test_s2       <- cbind(y=test_data_lapsed$early_lapse,  enc_s2$test)

make_foldid <- function(train_df, folds) {
  fid <- integer(nrow(train_df))
  for (k in seq_len(nrow(folds))) {
    va_ids <- unique(assessment(folds$splits[[k]])$id)
    fid[train_df$id %in% va_ids] <- k
  }
  fid
}
foldid_s1 <- make_foldid(train_data,        cv_folds)
foldid_cx <- make_foldid(train_data,        cv_folds)
foldid_s2 <- make_foldid(train_data_lapsed, cv_folds)

progress("Preprocessing complete.")

# ==============================================================================
# SECTION 5: FULL BASELINE MODELS
# Fit once on the full training set.
#   (a) Starting point for full-training stepAIC
#   (b) Baseline row in each comparison table
# ==============================================================================

progress("Fitting full baseline models...")

full_s1  <- glm(y ~ ., data = train_mod_lap, family = binomial())
full_cox  <- coxph(Surv(time, event) ~ ., data = train_mod_srv)
full_s2  <- glm(y ~ ., data = train_s2,    family = binomial())

known_s1  <- get_predictor_cols(train_data,        "lapse_flag")
known_cx  <- get_predictor_cols(train_data,        "days_to_event")
known_s2  <- get_predictor_cols(train_data_lapsed, "early_lapse")

progress("Full baseline models fitted.")

# ==============================================================================
# SECTION 6: CV EVALUATION HELPERS (used by baseline only)
# For Sections 8-10, selection runs inside each fold rather than evaluating
# a fixed feature set. cv_eval_binary / cv_eval_cox are retained for the
# baseline rows and can be reused wherever a fixed feature set is needed.
# ==============================================================================

cv_eval_binary <- function(cv_folds, train_data, vars, target,
                           lapsed_only = FALSE) {
  lapply(cv_folds$splits, function(sp) {
    tr <- analysis(sp); va <- assessment(sp)
    if (lapsed_only) {
      tr <- tr %>% filter(lapse_flag == 1) %>%
        mutate(early_lapse = as.integer(days_to_event <= LANDMARK_T1))
      va <- va %>% filter(lapse_flag == 1) %>%
        mutate(early_lapse = as.integer(days_to_event <= LANDMARK_T1))
      if (nrow(tr) < 10 || nrow(va) < 2 ||
          length(unique(tr[[target]])) < 2) return(NULL)
    }
    tryCatch({
      enc        <- encode_predictors(tr, va, target)
      avail      <- intersect(vars, names(enc$train))
      X_tr       <- cbind(y = tr[[target]], enc$train[, avail, drop = FALSE])
      X_va       <- enc$test[, avail, drop = FALSE]
      m          <- glm(y ~ ., data = X_tr, family = binomial())
      p          <- predict(m, X_va, type = "response")
      fm         <- binary_lapse_metrics(va[[target]], p)
      fm$AIC     <- AIC(m)
      fm
    }, error = function(e) NULL)
  })
}

cv_eval_cox <- function(cv_folds, train_data, vars) {
  lapply(cv_folds$splits, function(sp) {
    tr <- analysis(sp); va <- assessment(sp)
    tryCatch({
      enc    <- encode_predictors(tr, va, "days_to_event")
      avail  <- intersect(vars, names(enc$train))
      tr_mod <- cbind(time=tr$days_to_event, event=tr$lapse_flag,
                      enc$train[, avail, drop=FALSE])
      va_mod <- cbind(time=va$days_to_event, event=va$lapse_flag,
                      enc$test[,  avail, drop=FALSE])
      m  <- coxph(Surv(time, event) ~ ., data = tr_mod)
      lp <- predict(m, va_mod, type = "lp")
      sp_va <- get_surv_from_cox(m, lp)
      bm_W  <- binary_lapse_metrics(va$lapse_flag, 1 - sp_va$S_T2)
      bm_WN <- tryCatch({
        lv <- which(va_mod$event == 1L)
        if (length(lv) < 2) stop("too few")
        ev <- as.integer(va_mod$time[lv] <= LANDMARK_T1)
        if (length(unique(ev)) < 2) stop("no variation")
        cp <- pmin((1-sp_va$S_T1[lv])/pmax(1-sp_va$S_T2[lv],1e-10),1)
        binary_lapse_metrics(ev, cp)
      }, error = function(e)
        data.frame(R2=NA,RMSE=NA,MAE=NA,Brier=NA,Binom_Dev=NA,
                   AUC=NA,Gini=NA,Accuracy=NA,Precision=NA,
                   Recall=NA,Specificity=NA,F1=NA))
      data.frame(
        C_Index  = cindex_score(va_mod$time, va_mod$event, lp),
        tdAUC_T1 = td_auc_at(va_mod$time, va_mod$event, lp, LANDMARK_T1),
        tdAUC_T2 = td_auc_at(va_mod$time, va_mod$event, lp, LANDMARK_T2),
        Brier_T1 = ipcw_brier_score(tr_mod$time,tr_mod$event,
                                    va_mod$time,va_mod$event,sp_va$S_T1,LANDMARK_T1),
        Brier_T2 = ipcw_brier_score(tr_mod$time,tr_mod$event,
                                    va_mod$time,va_mod$event,sp_va$S_T2,LANDMARK_T2),
        IBS      = integrated_brier_score(tr_mod$time,tr_mod$event,
                                          va_mod$time,va_mod$event,m,lp),
        AIC      = AIC(m),
        W_AUC    = bm_W$AUC,   W_Brier  = bm_W$Brier,  W_Gini  = bm_W$Gini,
        W_F1     = bm_W$F1,
        WN_AUC   = bm_WN$AUC,  WN_Brier = bm_WN$Brier, WN_Gini = bm_WN$Gini,
        WN_F1    = bm_WN$F1
      )
    }, error = function(e) NULL)
  })
}

test_eval_binary <- function(model, test_data, target) {
  p <- predict(model, test_data, type = "response")
  binary_lapse_metrics(test_data[[target]], p)
}

make_row_binary <- function(method, n_vars, cv_sum, test_m) {
  data.frame(
    Method        = method,
    N_Vars        = n_vars,
    CV_AUC_Mean   = round(cv_sum["Mean","AUC"],   4),
    CV_AUC_SE     = round(cv_sum["SE",  "AUC"],   4),
    CV_Brier_Mean = round(cv_sum["Mean","Brier"],  4),
    CV_Brier_SE   = round(cv_sum["SE",  "Brier"],  4),
    CV_F1_Mean    = round(cv_sum["Mean","F1"],     4),
    CV_F1_SE      = round(cv_sum["SE",  "F1"],     4),
    Test_AUC      = round(test_m$AUC,              4),
    Test_Brier    = round(test_m$Brier,             4),
    Test_F1       = round(test_m$F1,               4),
    stringsAsFactors = FALSE
  )
}

make_row_cox <- function(method, n_vars, cv_sum, ci_te, ibs_te, wauc_te, wibs_te) {
  data.frame(
    Method           = method,
    N_Vars           = n_vars,
    CV_CIndex_Mean   = round(cv_sum["Mean","C_Index"], 4),
    CV_CIndex_SE     = round(cv_sum["SE",  "C_Index"], 4),
    CV_IBS_Mean      = round(cv_sum["Mean","IBS"],     4),
    CV_IBS_SE        = round(cv_sum["SE",  "IBS"],     4),
    CV_W_AUC_Mean    = round(cv_sum["Mean","W_AUC"],   4),
    CV_W_AUC_SE      = round(cv_sum["SE",  "W_AUC"],   4),
    CV_WN_AUC_Mean   = round(cv_sum["Mean","WN_AUC"],  4),
    CV_WN_AUC_SE     = round(cv_sum["SE",  "WN_AUC"],  4),
    Test_CIndex      = round(ci_te,                     4),
    Test_IBS         = round(ibs_te,                    4),
    Test_W_AUC       = round(wauc_te,                   4),
    Test_WN_AUC      = round(wibs_te,                   4),
    stringsAsFactors = FALSE
  )
}

# Per-fold Cox metric helper — reused across Boruta and Elastic Net Cox loops
fold_cox_metrics <- function(tr, va, tr_sel, va_sel) {
  m  <- tryCatch(coxph(Surv(time,event) ~ ., data=tr_sel), error=function(e) NULL)
  if (is.null(m)) return(NULL)
  lp    <- predict(m, va_sel, type="lp")
  sp_va <- get_surv_from_cox(m, lp)
  bm_W  <- binary_lapse_metrics(va$lapse_flag, 1 - sp_va$S_T2)
  bm_WN <- tryCatch({
    lv <- which(va_sel$event == 1L)
    if (length(lv) < 2) stop("too few")
    ev <- as.integer(va_sel$time[lv] <= LANDMARK_T1)
    if (length(unique(ev)) < 2) stop("no variation")
    cp <- pmin((1-sp_va$S_T1[lv])/pmax(1-sp_va$S_T2[lv],1e-10),1)
    binary_lapse_metrics(ev, cp)
  }, error = function(e)
    data.frame(R2=NA,RMSE=NA,MAE=NA,Brier=NA,Binom_Dev=NA,
               AUC=NA,Gini=NA,Accuracy=NA,Precision=NA,
               Recall=NA,Specificity=NA,F1=NA))
  data.frame(
    C_Index  = cindex_score(va_sel$time, va_sel$event, lp),
    tdAUC_T1 = td_auc_at(va_sel$time, va_sel$event, lp, LANDMARK_T1),
    tdAUC_T2 = td_auc_at(va_sel$time, va_sel$event, lp, LANDMARK_T2),
    Brier_T1 = ipcw_brier_score(tr_sel$time,tr_sel$event,
                                va_sel$time,va_sel$event,sp_va$S_T1,LANDMARK_T1),
    Brier_T2 = ipcw_brier_score(tr_sel$time,tr_sel$event,
                                va_sel$time,va_sel$event,sp_va$S_T2,LANDMARK_T2),
    IBS      = integrated_brier_score(tr_sel$time,tr_sel$event,
                                      va_sel$time,va_sel$event,m,lp),
    AIC      = AIC(m),
    W_AUC    = bm_W$AUC,   W_Brier  = bm_W$Brier,  W_Gini  = bm_W$Gini,
    W_F1     = bm_W$F1,
    WN_AUC   = bm_WN$AUC,  WN_Brier = bm_WN$Brier, WN_Gini = bm_WN$Gini,
    WN_F1    = bm_WN$F1
  )
}

# ==============================================================================
# SECTION 7: BASELINE CV EVALUATION
# ==============================================================================

progress("Evaluating baseline models via CV...")

vars_full_s1   <- names(enc_s1$train)
cv_raw_base_s1 <- cv_eval_binary(cv_folds, train_data, vars_full_s1, "lapse_flag")
cv_raw_base_s1 <- Filter(Negate(is.null), cv_raw_base_s1)
cv_sum_base_s1 <- cv_summary(do.call(rbind, cv_raw_base_s1))
test_base_s1   <- test_eval_binary(full_s1, test_mod_lap, "y")
row_base_s1    <- make_row_binary("Baseline (Full)", ncol(enc_s1$train),
                                  cv_sum_base_s1, test_base_s1)
progress("  Stage 1 baseline CV done.")

vars_full_cx   <- names(enc_cx$train)
cv_raw_base_cx <- cv_eval_cox(cv_folds, train_data, vars_full_cx)
cv_raw_base_cx <- Filter(Negate(is.null), cv_raw_base_cx)
cv_sum_base_cx <- cv_summary(do.call(rbind, cv_raw_base_cx))
lp_test_cx     <- predict(full_cox, test_mod_srv, type = "lp")
sp_test_cx     <- get_surv_from_cox(full_cox, lp_test_cx)
ci_test_cx     <- cindex_score(test_mod_srv$time, test_mod_srv$event, lp_test_cx)
ibs_test_cx    <- integrated_brier_score(train_mod_srv$time, train_mod_srv$event,
                                         test_mod_srv$time,  test_mod_srv$event,
                                         full_cox, lp_test_cx)
wauc_test_cx   <- auc_score(test_data$lapse_flag, 1 - sp_test_cx$S_T2)
lapsed_te_idx  <- which(test_mod_srv$event == 1L)
cox_cond_te    <- pmin((1-sp_test_cx$S_T1[lapsed_te_idx]) /
                         pmax(1-sp_test_cx$S_T2[lapsed_te_idx],1e-10), 1)
early_te_obs   <- as.integer(test_mod_srv$time[lapsed_te_idx] <= LANDMARK_T1)
wnauc_test_cx  <- auc_score(early_te_obs, cox_cond_te)
row_base_cx    <- make_row_cox("Baseline (Full)", ncol(enc_cx$train),
                               cv_sum_base_cx, ci_test_cx, ibs_test_cx,
                               wauc_test_cx, wnauc_test_cx)
progress("  Cox baseline CV done.")

vars_full_s2   <- names(enc_s2$train)
cv_raw_base_s2 <- cv_eval_binary(cv_folds, train_data_lapsed, vars_full_s2,
                                 "early_lapse", lapsed_only = TRUE)
cv_raw_base_s2 <- Filter(Negate(is.null), cv_raw_base_s2)
cv_sum_base_s2 <- cv_summary(do.call(rbind, cv_raw_base_s2))
test_base_s2   <- test_eval_binary(full_s2, test_s2, "y")
row_base_s2    <- make_row_binary("Baseline (Full)", ncol(enc_s2$train),
                                  cv_sum_base_s2, test_base_s2)
progress("  Stage 2 baseline CV done.")

save_results(row_base_s1, "baseline", "stage1")
save_results(row_base_cx, "baseline", "cox")
save_results(row_base_s2, "baseline", "stage2")
progress("Baseline CV evaluation complete.")

# ==============================================================================
# SECTION 8: STEPWISE AIC
# stepAIC runs INSIDE each CV fold (fold-blind selection).
# Full-training stepAIC produces the final feature set for test evaluation.
# ==============================================================================

progress("=== STEPWISE AIC ===")

# ------------------------------------------------------------------------------
# 8.1 Stage 1 Stepwise — binary logistic
# ------------------------------------------------------------------------------
progress("  Stage 1: running per-fold stepAIC CV...")

cv_step_s1_raw <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  cat(sprintf("    [Stepwise S1] Fold %d/5 ... ", i)); flush.console()
  tr <- analysis(sp); va <- assessment(sp)
  tryCatch({
    enc    <- encode_predictors(tr, va, "lapse_flag")
    full_m <- glm(y ~ ., data=cbind(y=tr$lapse_flag, enc$train), family=binomial())
    step_m <- stepAIC(full_m, direction="both", trace=0)
    sel    <- intersect(attr(terms(step_m),"term.labels"),
                        intersect(names(enc$train), names(enc$test)))
    if (length(sel) == 0) { cat("skipped (no features)\n"); return(NULL) }
    m   <- glm(y ~ ., data=cbind(y=tr$lapse_flag, enc$train[,sel,drop=FALSE]),
               family=binomial())
    p   <- predict(m, enc$test[,sel,drop=FALSE], type="response")
    fm  <- binary_lapse_metrics(va$lapse_flag, p)
    fm$AIC <- AIC(m)
    cat(sprintf("done  vars=%d  AUC=%.4f\n", length(sel), fm$AUC)); flush.console()
    fm
  }, error=function(e) { cat("error\n"); NULL })
})
cv_raw_step_s1 <- Filter(Negate(is.null), cv_step_s1_raw)
cv_sum_step_s1 <- cv_summary(do.call(rbind, cv_raw_step_s1))

# Full-training stepAIC → final feature set
progress("  Stage 1: running stepAIC on full training set...")
step_s1            <- stepAIC(full_s1, direction="both", trace=0)
step_vars_s1_dummy <- attr(terms(step_s1), "term.labels")
step_vars_s1_orig  <- unique(strip_suffix(step_vars_s1_dummy, known_s1))
progress(sprintf("  Stage 1: %d original variables selected.", length(step_vars_s1_orig)))
save_results(data.frame(variable=step_vars_s1_orig, stringsAsFactors=FALSE),
             "stepwise", "stage1", "selected_vars")
test_step_s1 <- test_eval_binary(step_s1, test_mod_lap, "y")
row_step_s1  <- make_row_binary("Stepwise AIC", length(step_vars_s1_orig),
                                cv_sum_step_s1, test_step_s1)
save_results(do.call(rbind, cv_raw_step_s1), "stepwise", "stage1", "cv_results")
save_results(row_step_s1,                    "stepwise", "stage1", "metrics")
progress("  Stage 1 Stepwise complete.")

# ------------------------------------------------------------------------------
# 8.2 Cox Stepwise
# ------------------------------------------------------------------------------
progress("  Cox: running per-fold stepAIC CV...")

cv_step_cx_raw <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  cat(sprintf("    [Stepwise Cox] Fold %d/5 ... ", i)); flush.console()
  tr <- analysis(sp); va <- assessment(sp)
  tryCatch({
    enc    <- encode_predictors(tr, va, "days_to_event")
    tr_mod <- cbind(time=tr$days_to_event, event=tr$lapse_flag, enc$train)
    full_m <- coxph(Surv(time,event) ~ ., data=tr_mod)
    step_m <- stepAIC(full_m, direction="both", trace=0)
    sel    <- intersect(attr(terms(step_m),"term.labels"),
                        intersect(names(enc$train), names(enc$test)))
    if (length(sel) == 0) { cat("skipped (no features)\n"); return(NULL) }
    tr_sel <- cbind(time=tr$days_to_event, event=tr$lapse_flag,
                    enc$train[,sel,drop=FALSE])
    va_sel <- cbind(time=va$days_to_event, event=va$lapse_flag,
                    enc$test[,sel,drop=FALSE])
    res <- fold_cox_metrics(tr, va, tr_sel, va_sel)
    if (!is.null(res))
      cat(sprintf("done  vars=%d  CIndex=%.4f\n", length(sel), res$C_Index))
    else cat("skipped\n")
    flush.console()
    res
  }, error=function(e) { cat("error\n"); NULL })
})
cv_raw_step_cx <- Filter(Negate(is.null), cv_step_cx_raw)
cv_sum_step_cx <- cv_summary(do.call(rbind, cv_raw_step_cx))

progress("  Cox: running stepAIC on full training set...")
step_cox           <- stepAIC(full_cox, direction="both", trace=0)
step_vars_cx_dummy <- attr(terms(step_cox), "term.labels")
step_vars_cx_orig  <- unique(strip_suffix(step_vars_cx_dummy, known_cx))
progress(sprintf("  Cox: %d original variables selected.", length(step_vars_cx_orig)))
save_results(data.frame(variable=step_vars_cx_orig, stringsAsFactors=FALSE),
             "stepwise", "cox", "selected_vars")
lp_step_cx_te  <- predict(step_cox, test_mod_srv, type="lp")
sp_step_cx_te  <- get_surv_from_cox(step_cox, lp_step_cx_te)
ci_step_cx     <- cindex_score(test_mod_srv$time, test_mod_srv$event, lp_step_cx_te)
ibs_step_cx    <- integrated_brier_score(train_mod_srv$time, train_mod_srv$event,
                                         test_mod_srv$time,  test_mod_srv$event,
                                         step_cox, lp_step_cx_te)
wauc_step_cx   <- auc_score(test_data$lapse_flag, 1-sp_step_cx_te$S_T2)
cond_step_cx   <- pmin((1-sp_step_cx_te$S_T1[lapsed_te_idx]) /
                         pmax(1-sp_step_cx_te$S_T2[lapsed_te_idx],1e-10), 1)
wnauc_step_cx  <- auc_score(early_te_obs, cond_step_cx)
row_step_cx    <- make_row_cox("Stepwise AIC", length(step_vars_cx_orig),
                               cv_sum_step_cx, ci_step_cx, ibs_step_cx,
                               wauc_step_cx, wnauc_step_cx)
save_results(do.call(rbind, cv_raw_step_cx), "stepwise", "cox", "cv_results")
save_results(row_step_cx,                    "stepwise", "cox", "metrics")
progress("  Cox Stepwise complete.")

# ------------------------------------------------------------------------------
# 8.3 Stage 2 Stepwise — binary logistic, lapsed only
# ------------------------------------------------------------------------------
progress("  Stage 2: running per-fold stepAIC CV...")

cv_step_s2_raw <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  cat(sprintf("    [Stepwise S2] Fold %d/5 ... ", i)); flush.console()
  tr <- analysis(sp) %>% filter(lapse_flag==1) %>%
        mutate(early_lapse=as.integer(days_to_event<=LANDMARK_T1))
  va <- assessment(sp) %>% filter(lapse_flag==1) %>%
        mutate(early_lapse=as.integer(days_to_event<=LANDMARK_T1))
  if (nrow(tr)<10 || nrow(va)<2 || length(unique(tr$early_lapse))<2) {
    cat("skipped (insufficient lapsed rows)\n"); return(NULL)
  }
  tryCatch({
    enc    <- encode_predictors(tr, va, "early_lapse")
    full_m <- glm(y ~ ., data=cbind(y=tr$early_lapse, enc$train), family=binomial())
    step_m <- stepAIC(full_m, direction="both", trace=0)
    sel    <- intersect(attr(terms(step_m),"term.labels"),
                        intersect(names(enc$train), names(enc$test)))
    if (length(sel)==0) { cat("skipped (no features)\n"); return(NULL) }
    m  <- glm(y ~ ., data=cbind(y=tr$early_lapse, enc$train[,sel,drop=FALSE]),
              family=binomial())
    p  <- predict(m, enc$test[,sel,drop=FALSE], type="response")
    fm <- binary_lapse_metrics(va$early_lapse, p)
    fm$AIC <- AIC(m)
    cat(sprintf("done  vars=%d  AUC=%.4f\n", length(sel), fm$AUC)); flush.console()
    fm
  }, error=function(e) { cat("error\n"); NULL })
})
cv_raw_step_s2 <- Filter(Negate(is.null), cv_step_s2_raw)
cv_sum_step_s2 <- cv_summary(do.call(rbind, cv_raw_step_s2))

progress("  Stage 2: running stepAIC on full training set...")
step_s2            <- stepAIC(full_s2, direction="both", trace=0)
step_vars_s2_dummy <- attr(terms(step_s2), "term.labels")
step_vars_s2_orig  <- unique(strip_suffix(step_vars_s2_dummy, known_s2))
progress(sprintf("  Stage 2: %d original variables selected.", length(step_vars_s2_orig)))
save_results(data.frame(variable=step_vars_s2_orig, stringsAsFactors=FALSE),
             "stepwise", "stage2", "selected_vars")
test_step_s2 <- test_eval_binary(step_s2, test_s2, "y")
row_step_s2  <- make_row_binary("Stepwise AIC", length(step_vars_s2_orig),
                                cv_sum_step_s2, test_step_s2)
save_results(do.call(rbind, cv_raw_step_s2), "stepwise", "stage2", "cv_results")
save_results(row_step_s2,                    "stepwise", "stage2", "metrics")
progress("  Stage 2 Stepwise complete.")

# ==============================================================================
# SECTION 9: BORUTA
# Boruta runs INSIDE each CV fold (fold-blind selection).
# doTrace=2  : numbered per-run lines + summary blocks (native output, no wrapper)
# verbose=FALSE : suppresses ranger's own per-tree progress
# Full-training Boruta produces the final feature set for test evaluation.
# BORUTA_RUNS = 200 (kept from original; Cox model needs more iterations)
# ==============================================================================

progress("=== BORUTA ===")
BORUTA_RUNS <- 200

# ------------------------------------------------------------------------------
# 9.1 Stage 1 Boruta — classification RF on factor(lapse_flag)
# ------------------------------------------------------------------------------
progress(sprintf("  Stage 1: running per-fold Boruta CV (maxRuns=%d)...", BORUTA_RUNS))

cv_bor_s1_raw <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  progress(sprintf("\n  \u2500\u2500 [Boruta S1] Fold %d/5 \u2500\u2500", i))
  t_fold <- proc.time()[["elapsed"]]
  tr <- analysis(sp); va <- assessment(sp)
  tryCatch({
    enc <- encode_predictors(tr, va, "lapse_flag")
    set.seed(42 + i)
    bor <- Boruta(
      x       = enc$train,
      y       = factor(tr$lapse_flag),
      maxRuns = BORUTA_RUNS,
      doTrace = 2,
      verbose = FALSE
    )
    bor       <- TentativeRoughFix(bor)
    confirmed <- getSelectedAttributes(bor, withTentative=FALSE)
    if (length(confirmed)==0) confirmed <- names(enc$train)
    sel <- intersect(confirmed, intersect(names(enc$train), names(enc$test)))
    if (length(sel)==0) return(NULL)
    m  <- glm(y ~ ., data=cbind(y=tr$lapse_flag, enc$train[,sel,drop=FALSE]),
              family=binomial())
    p  <- predict(m, enc$test[,sel,drop=FALSE], type="response")
    fm <- binary_lapse_metrics(va$lapse_flag, p)
    fm$AIC <- AIC(m)
    elapsed <- proc.time()[["elapsed"]] - t_fold
    progress(sprintf("  -> Fold %d/5 done (%.1f mins) \u2014 %d vars | AUC=%.4f",
                     i, elapsed/60, length(sel), fm$AUC))
    fm
  }, error=function(e) { progress(sprintf("  Fold %d error: %s", i, conditionMessage(e))); NULL })
})
cv_raw_bor_s1 <- Filter(Negate(is.null), cv_bor_s1_raw)
cv_sum_bor_s1 <- cv_summary(do.call(rbind, cv_raw_bor_s1))

# Full-training Boruta
progress(sprintf("\n  Stage 1: running Boruta on full training set (maxRuns=%d)...",
                 BORUTA_RUNS))
cat("Progress (doTrace=2):\n\n")
set.seed(42)
boruta_s1 <- Boruta(
  x       = as.data.frame(enc_s1$train),
  y       = factor(train_data$lapse_flag),
  maxRuns = BORUTA_RUNS,
  doTrace = 2,
  verbose = FALSE
)
boruta_s1_final <- TentativeRoughFix(boruta_s1)
confirmed_s1    <- names(boruta_s1_final$finalDecision[
                     boruta_s1_final$finalDecision == "Confirmed"])
orig_s1         <- unique(strip_suffix(confirmed_s1, known_s1))
progress(sprintf("  Stage 1: %d confirmed variables.", length(orig_s1)))
save_results(data.frame(variable=orig_s1, stringsAsFactors=FALSE),
             "boruta", "stage1", "selected_vars")
saveRDS(boruta_s1_final, here("fs_boruta_stage1_boruta_object.rds"))

bor_s1_refit <- glm(
  as.formula(paste("y ~", paste(confirmed_s1, collapse=" + "))),
  data=train_mod_lap[, c("y", confirmed_s1)], family=binomial()
)
test_bor_s1 <- test_eval_binary(bor_s1_refit, test_mod_lap, "y")
row_bor_s1  <- make_row_binary("Boruta", length(orig_s1), cv_sum_bor_s1, test_bor_s1)
save_results(do.call(rbind, cv_raw_bor_s1), "boruta", "stage1", "cv_results")
save_results(row_bor_s1,                    "boruta", "stage1", "metrics")
progress("  Stage 1 Boruta complete.")

# ------------------------------------------------------------------------------
# 9.2 Cox Boruta — martingale residuals as continuous surrogate response
# ------------------------------------------------------------------------------
progress(sprintf("  Cox: running per-fold Boruta CV (maxRuns=%d)...", BORUTA_RUNS))

cv_bor_cx_raw <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  progress(sprintf("\n  \u2500\u2500 [Boruta Cox] Fold %d/5 \u2500\u2500", i))
  t_fold <- proc.time()[["elapsed"]]
  tr <- analysis(sp); va <- assessment(sp)
  tryCatch({
    enc    <- encode_predictors(tr, va, "days_to_event")
    tr_mod <- cbind(time=tr$days_to_event, event=tr$lapse_flag, enc$train)
    # Null Cox → martingale residuals as Boruta target
    null_m <- coxph(Surv(time,event) ~ 1, data=tr_mod)
    mart   <- residuals(null_m, type="martingale")
    set.seed(42 + i)
    bor <- Boruta(
      x       = enc$train,
      y       = mart,
      maxRuns = BORUTA_RUNS,
      doTrace = 2,
      verbose = FALSE
    )
    bor       <- TentativeRoughFix(bor)
    confirmed <- getSelectedAttributes(bor, withTentative=FALSE)
    if (length(confirmed)==0) confirmed <- names(enc$train)
    sel <- intersect(confirmed, intersect(names(enc$train), names(enc$test)))
    if (length(sel)==0) return(NULL)
    tr_sel <- cbind(time=tr$days_to_event, event=tr$lapse_flag,
                    enc$train[,sel,drop=FALSE])
    va_sel <- cbind(time=va$days_to_event, event=va$lapse_flag,
                    enc$test[,sel,drop=FALSE])
    res <- fold_cox_metrics(tr, va, tr_sel, va_sel)
    elapsed <- proc.time()[["elapsed"]] - t_fold
    if (!is.null(res))
      progress(sprintf("  -> Fold %d/5 done (%.1f mins) \u2014 %d vars | C-Index=%.4f",
                       i, elapsed/60, length(sel), res$C_Index))
    res
  }, error=function(e) { progress(sprintf("  Fold %d error: %s", i, conditionMessage(e))); NULL })
})
cv_raw_bor_cx <- Filter(Negate(is.null), cv_bor_cx_raw)
cv_sum_bor_cx <- cv_summary(do.call(rbind, cv_raw_bor_cx))

# Full-training Cox Boruta
progress(sprintf("\n  Cox: computing martingale residuals and running Boruta (maxRuns=%d)...",
                 BORUTA_RUNS))
cat("Progress (doTrace=2):\n\n")
null_cox   <- coxph(Surv(time,event) ~ 1, data=train_mod_srv)
mart_resid <- residuals(null_cox, type="martingale")
set.seed(42)
boruta_cx <- Boruta(
  x       = as.data.frame(enc_cx$train),
  y       = mart_resid,
  maxRuns = BORUTA_RUNS,
  doTrace = 2,
  verbose = FALSE
)
boruta_cx_final <- TentativeRoughFix(boruta_cx)
confirmed_cx    <- names(boruta_cx_final$finalDecision[
                     boruta_cx_final$finalDecision == "Confirmed"])
orig_cx         <- unique(strip_suffix(confirmed_cx, known_cx))
progress(sprintf("  Cox: %d confirmed variables.", length(orig_cx)))
save_results(data.frame(variable=orig_cx, stringsAsFactors=FALSE),
             "boruta", "cox", "selected_vars")
saveRDS(boruta_cx_final, here("fs_boruta_cox_boruta_object.rds"))

bor_cx_refit <- coxph(
  as.formula(paste("Surv(time,event) ~", paste(confirmed_cx, collapse=" + "))),
  data=train_mod_srv[, c("time","event", confirmed_cx)]
)
lp_bor_cx_te  <- predict(bor_cx_refit, test_mod_srv, type="lp")
sp_bor_cx_te  <- get_surv_from_cox(bor_cx_refit, lp_bor_cx_te)
ci_bor_cx     <- cindex_score(test_mod_srv$time, test_mod_srv$event, lp_bor_cx_te)
ibs_bor_cx    <- integrated_brier_score(train_mod_srv$time, train_mod_srv$event,
                                        test_mod_srv$time,  test_mod_srv$event,
                                        bor_cx_refit, lp_bor_cx_te)
wauc_bor_cx   <- auc_score(test_data$lapse_flag, 1-sp_bor_cx_te$S_T2)
cond_bor_cx   <- pmin((1-sp_bor_cx_te$S_T1[lapsed_te_idx]) /
                        pmax(1-sp_bor_cx_te$S_T2[lapsed_te_idx],1e-10), 1)
wnauc_bor_cx  <- auc_score(early_te_obs, cond_bor_cx)
row_bor_cx    <- make_row_cox("Boruta", length(orig_cx), cv_sum_bor_cx,
                              ci_bor_cx, ibs_bor_cx, wauc_bor_cx, wnauc_bor_cx)
save_results(do.call(rbind, cv_raw_bor_cx), "boruta", "cox", "cv_results")
save_results(row_bor_cx,                    "boruta", "cox", "metrics")
progress("  Cox Boruta complete.")

# ------------------------------------------------------------------------------
# 9.3 Stage 2 Boruta — classification RF, lapsed-only subset
# ------------------------------------------------------------------------------
progress(sprintf("  Stage 2: running per-fold Boruta CV (maxRuns=%d)...", BORUTA_RUNS))

cv_bor_s2_raw <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  progress(sprintf("\n  \u2500\u2500 [Boruta S2] Fold %d/5 \u2500\u2500", i))
  t_fold <- proc.time()[["elapsed"]]
  tr <- analysis(sp) %>% filter(lapse_flag==1) %>%
        mutate(early_lapse=as.integer(days_to_event<=LANDMARK_T1))
  va <- assessment(sp) %>% filter(lapse_flag==1) %>%
        mutate(early_lapse=as.integer(days_to_event<=LANDMARK_T1))
  if (nrow(tr)<10 || nrow(va)<2 || length(unique(tr$early_lapse))<2) {
    progress(sprintf("  Fold %d skipped (insufficient lapsed rows)", i)); return(NULL)
  }
  tryCatch({
    enc <- encode_predictors(tr, va, "early_lapse")
    set.seed(42 + i)
    bor <- Boruta(
      x       = enc$train,
      y       = factor(tr$early_lapse),
      maxRuns = BORUTA_RUNS,
      doTrace = 2,
      verbose = FALSE
    )
    bor       <- TentativeRoughFix(bor)
    confirmed <- getSelectedAttributes(bor, withTentative=FALSE)
    if (length(confirmed)==0) confirmed <- names(enc$train)
    sel <- intersect(confirmed, intersect(names(enc$train), names(enc$test)))
    if (length(sel)==0) return(NULL)
    m  <- glm(y ~ ., data=cbind(y=tr$early_lapse, enc$train[,sel,drop=FALSE]),
              family=binomial())
    p  <- predict(m, enc$test[,sel,drop=FALSE], type="response")
    fm <- binary_lapse_metrics(va$early_lapse, p)
    fm$AIC <- AIC(m)
    elapsed <- proc.time()[["elapsed"]] - t_fold
    progress(sprintf("  -> Fold %d/5 done (%.1f mins) \u2014 %d vars | AUC=%.4f",
                     i, elapsed/60, length(sel), fm$AUC))
    fm
  }, error=function(e) { progress(sprintf("  Fold %d error: %s", i, conditionMessage(e))); NULL })
})
cv_raw_bor_s2 <- Filter(Negate(is.null), cv_bor_s2_raw)

# Full-training Stage 2 Boruta
progress(sprintf("\n  Stage 2: running Boruta on full training set (maxRuns=%d)...",
                 BORUTA_RUNS))
cat("Progress (doTrace=2):\n\n")
set.seed(42)
boruta_s2 <- Boruta(
  x       = as.data.frame(enc_s2$train),
  y       = factor(train_data_lapsed$early_lapse),
  maxRuns = BORUTA_RUNS,
  doTrace = 2,
  verbose = FALSE
)
boruta_s2_final <- TentativeRoughFix(boruta_s2)
confirmed_s2    <- names(boruta_s2_final$finalDecision[
                     boruta_s2_final$finalDecision == "Confirmed"])
orig_s2         <- unique(strip_suffix(confirmed_s2, known_s2))
progress(sprintf("  Stage 2: %d confirmed variables.", length(orig_s2)))
save_results(data.frame(variable=orig_s2, stringsAsFactors=FALSE),
             "boruta", "stage2", "selected_vars")
saveRDS(boruta_s2_final, here("fs_boruta_stage2_boruta_object.rds"))

if (length(confirmed_s2) == 0) {
  progress("  WARNING: Zero variables confirmed for Stage 2.")
  row_bor_s2 <- data.frame(Method="Boruta", N_Vars=0,
                           CV_AUC_Mean=NA, CV_AUC_SE=NA,
                           CV_Brier_Mean=NA, CV_Brier_SE=NA,
                           CV_F1_Mean=NA, CV_F1_SE=NA,
                           Test_AUC=NA, Test_Brier=NA, Test_F1=NA,
                           stringsAsFactors=FALSE)
  save_results(row_bor_s2, "boruta", "stage2", "metrics")
} else if (length(cv_raw_bor_s2) == 0) {
  progress("  WARNING: All Stage 2 Boruta CV folds skipped — no CV summary available.")
  row_bor_s2 <- data.frame(Method="Boruta", N_Vars=length(orig_s2),
                           CV_AUC_Mean=NA, CV_AUC_SE=NA,
                           CV_Brier_Mean=NA, CV_Brier_SE=NA,
                           CV_F1_Mean=NA, CV_F1_SE=NA,
                           Test_AUC=NA, Test_Brier=NA, Test_F1=NA,
                           stringsAsFactors=FALSE)
  save_results(row_bor_s2, "boruta", "stage2", "metrics")
} else {
  cv_sum_bor_s2 <- cv_summary(do.call(rbind, cv_raw_bor_s2))
  safe_s2            <- paste0("`", confirmed_s2, "`")
  bor_s2_refit       <- glm(
    as.formula(paste("y ~", paste(safe_s2, collapse=" + "))),
    data=train_s2[, c("y", confirmed_s2)], family=binomial()
  )
  test_bor_s2  <- test_eval_binary(bor_s2_refit, test_s2, "y")
  row_bor_s2   <- make_row_binary("Boruta", length(orig_s2), cv_sum_bor_s2, test_bor_s2)
  save_results(do.call(rbind, cv_raw_bor_s2), "boruta", "stage2", "cv_results")
  save_results(row_bor_s2,                    "boruta", "stage2", "metrics")
}
progress("  Stage 2 Boruta complete.")

# ==============================================================================
# SECTION 10: ELASTIC NET
# Architecture (fold-blind — consistent with frequency/severity files):
#   Step 1 — Outer fold × ALPHA_GRID: cv.glmnet inside each fold (inner 5-fold
#             for lambda), predict on fold validation, record primary metric.
#             Best alpha = maximum mean outer-fold AUC (binary) / C-Index (Cox).
#   Step 2 — Per-fold CV with best_alpha: cv.glmnet inside each fold with
#             best_alpha → lambda.min and lambda.1se variable sets → unpenalised
#             GLM/Cox → full metric suite for comparison table.
#   Step 3 — Full-training: cv.glmnet with best_alpha + foldid → lambda.min
#             and lambda.1se → unpenalised refit → test evaluation.
# ==============================================================================

progress("=== ELASTIC NET ===")
ALPHA_GRID <- seq(0, 1, by = 0.1)   # 11 values from Ridge to LASSO

extract_nonzero_vars <- function(cv_fit, s = "lambda.min") {
  cf  <- coef(cv_fit, s = s)
  nms <- rownames(cf)[cf[, 1] != 0 & rownames(cf) != "(Intercept)"]
  nms
}
unpen_glm <- function(train_df, vars, target) {
  avail <- intersect(vars, names(train_df))
  glm(as.formula(paste(target, "~", paste(avail, collapse=" + "))),
      data=train_df[, c(target, avail)], family=binomial())
}
unpen_cox <- function(train_mod, vars) {
  avail <- intersect(vars, names(train_mod))
  coxph(as.formula(paste("Surv(time,event) ~", paste(avail, collapse=" + "))),
        data=train_mod[, c("time","event", avail)])
}

# ==============================================================================
# 10.1 Stage 1 Elastic Net
# ==============================================================================
progress("  Stage 1: outer fold x alpha grid (AUC criterion)...")

en_s1_outer <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  cat(sprintf("    [EN S1 alpha search] Fold %d/5 ... ", i)); flush.console()
  tr <- analysis(sp); va <- assessment(sp)
  tryCatch({
    enc     <- encode_predictors(tr, va, "lapse_flag")
    X_tr_f  <- as.matrix(enc$train); X_va_f <- as.matrix(enc$test)
    set.seed(42 + i)
    inner_clients     <- unique(tr$id)
    client_inner_fold <- setNames(
      sample(rep(1:5, length.out = length(inner_clients))),
      as.character(inner_clients)
    )
    inner_foldid <- client_inner_fold[as.character(tr$id)]
    
    fold_res <- lapply(ALPHA_GRID, function(a) {
      set.seed(42)
      cv_fit <- tryCatch(
        cv.glmnet(X_tr_f, tr$lapse_flag, family="binomial", alpha=a,
                  standardize=TRUE, foldid=inner_foldid, type.measure="auc"),
        error=function(e) NULL
      )
      if (is.null(cv_fit)) return(NULL)
      p_1se <- as.vector(predict(cv_fit, newx=X_va_f, s="lambda.1se", type="response"))
      data.frame(alpha=a, AUC_1se=auc_score(va$lapse_flag, p_1se))
    })
    res      <- do.call(rbind, Filter(Negate(is.null), fold_res))
    best_row <- res[which.max(res$AUC_1se), ]
    cat(sprintf("done  best alpha=%.1f  AUC=%.4f\n", best_row$alpha, best_row$AUC_1se))
    flush.console()
    res
  }, error=function(e) { cat("error\n"); NULL })
})
en_s1_outer_df <- do.call(rbind, Filter(Negate(is.null), en_s1_outer))
en_s1_alpha_sum <- aggregate(AUC_1se ~ alpha, data=en_s1_outer_df, FUN=mean)
best_alpha_s1   <- en_s1_alpha_sum$alpha[which.max(en_s1_alpha_sum$AUC_1se)]
progress(sprintf("  Stage 1: best alpha=%.1f", best_alpha_s1))

# Per-fold CV with best_alpha_s1 → lambda.min and lambda.1se metrics
progress(sprintf("  Stage 1: per-fold CV with best alpha=%.1f...", best_alpha_s1))

cv_en_s1_raw <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  cat(sprintf("    [EN S1 CV] Fold %d/5 ... ", i)); flush.console()
  tr <- analysis(sp); va <- assessment(sp)
  tryCatch({
    enc    <- encode_predictors(tr, va, "lapse_flag")
    X_tr_f <- as.matrix(enc$train); X_va_f <- as.matrix(enc$test)
    set.seed(42 + i)
    inner_clients     <- unique(tr$id)
    client_inner_fold <- setNames(
      sample(rep(1:5, length.out = length(inner_clients))),
      as.character(inner_clients)
    )
    inner_foldid <- client_inner_fold[as.character(tr$id)]
    set.seed(42)
    cv_fit <- cv.glmnet(X_tr_f, tr$lapse_flag, family="binomial",
                        alpha=best_alpha_s1, standardize=TRUE, foldid=inner_foldid, type.measure="auc")
    make_fold_row <- function(s_val) {
      sel <- intersect(
        extract_nonzero_vars(cv_fit, s_val),
        intersect(names(enc$train), names(enc$test))
      )
      if (length(sel)==0) return(NULL)
      m  <- glm(y~., data=cbind(y=tr$lapse_flag, enc$train[,sel,drop=FALSE]),
                family=binomial())
      p  <- predict(m, enc$test[,sel,drop=FALSE], type="response")
      fm <- binary_lapse_metrics(va$lapse_flag, p); fm$AIC <- AIC(m); fm
    }
    list(min=make_fold_row("lambda.min"), se1=make_fold_row("lambda.1se"))
  }, error=function(e) { cat("error\n"); list(min=NULL, se1=NULL) })
})
cv_raw_en_s1_min <- Filter(Negate(is.null), lapply(cv_en_s1_raw, `[[`, "min"))
cv_raw_en_s1_1se <- Filter(Negate(is.null), lapply(cv_en_s1_raw, `[[`, "se1"))
cat(sprintf("    done  %d/%d folds succeeded\n",
            length(cv_raw_en_s1_min), length(cv_folds$splits))); flush.console()
cv_sum_en_s1_min <- if (length(cv_raw_en_s1_min)>0) cv_summary(do.call(rbind, cv_raw_en_s1_min)) else NULL
cv_sum_en_s1_1se <- if (length(cv_raw_en_s1_1se)>0) cv_summary(do.call(rbind, cv_raw_en_s1_1se)) else NULL

# Full-training with best_alpha_s1
X_s1 <- as.matrix(enc_s1$train); y_s1 <- train_data$lapse_flag
set.seed(42)
best_en_s1 <- cv.glmnet(X_s1, y_s1, family="binomial", alpha=best_alpha_s1,
                         foldid=foldid_s1, type.measure="auc", standardize=TRUE)
vars_en_s1_min <- extract_nonzero_vars(best_en_s1, "lambda.min")
vars_en_s1_1se <- extract_nonzero_vars(best_en_s1, "lambda.1se")
orig_en_s1_min <- unique(strip_suffix(vars_en_s1_min, known_s1))
orig_en_s1_1se <- unique(strip_suffix(vars_en_s1_1se, known_s1))
en_s1_min_fit  <- unpen_glm(as.data.frame(train_mod_lap), vars_en_s1_min, "y")
en_s1_1se_fit  <- unpen_glm(as.data.frame(train_mod_lap), vars_en_s1_1se, "y")
test_en_s1_min <- test_eval_binary(en_s1_min_fit, test_mod_lap, "y")
test_en_s1_1se <- test_eval_binary(en_s1_1se_fit, test_mod_lap, "y")
row_en_s1_min  <- make_row_binary(sprintf("EN a=%.1f lambda.min", best_alpha_s1),
                                  length(orig_en_s1_min), cv_sum_en_s1_min, test_en_s1_min)
row_en_s1_1se  <- make_row_binary(sprintf("EN a=%.1f lambda.1se", best_alpha_s1),
                                  length(orig_en_s1_1se), cv_sum_en_s1_1se, test_en_s1_1se)
save_results(data.frame(variable=orig_en_s1_min), "elastic_min", "stage1", "selected_vars")
save_results(do.call(rbind,cv_raw_en_s1_min),     "elastic_min", "stage1", "cv_results")
save_results(row_en_s1_min,                        "elastic_min", "stage1", "metrics")
save_results(data.frame(variable=orig_en_s1_1se), "elastic_1se", "stage1", "selected_vars")
save_results(do.call(rbind,cv_raw_en_s1_1se),     "elastic_1se", "stage1", "cv_results")
save_results(row_en_s1_1se,                        "elastic_1se", "stage1", "metrics")
progress(sprintf("  Stage 1 EN complete: lambda.min=%d vars | lambda.1se=%d vars.",
                 length(orig_en_s1_min), length(orig_en_s1_1se)))

# ==============================================================================
# 10.2 Cox Elastic Net
# ==============================================================================
progress("  Cox: outer fold x alpha grid (C-Index criterion)...")

en_cx_outer <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  cat(sprintf("    [EN Cox alpha search] Fold %d/5 ... ", i)); flush.console()
  tr <- analysis(sp); va <- assessment(sp)
  tryCatch({
    enc    <- encode_predictors(tr, va, "days_to_event")
    X_tr_f <- as.matrix(enc$train); X_va_f <- as.matrix(enc$test)
    y_cx_f <- Surv(tr$days_to_event, tr$lapse_flag)
    set.seed(42 + i)
    inner_clients     <- unique(tr$id)
    client_inner_fold <- setNames(
      sample(rep(1:5, length.out = length(inner_clients))),
      as.character(inner_clients)
    )
    inner_foldid <- client_inner_fold[as.character(tr$id)]
    
    fold_res <- lapply(ALPHA_GRID, function(a) {
      set.seed(42)
      cv_fit <- tryCatch(
        cv.glmnet(X_tr_f, y_cx_f, family="cox", alpha=a,
                  standardize=TRUE, foldid=inner_foldid, type.measure="C"),
        error=function(e) NULL
      )
      if (is.null(cv_fit)) return(NULL)
      lp_1se <- as.vector(predict(cv_fit, newx=X_va_f, s="lambda.1se", type="link"))
      data.frame(alpha=a,
                 CIndex_1se=cindex_score(va$days_to_event, va$lapse_flag, lp_1se))
    })
    res      <- do.call(rbind, Filter(Negate(is.null), fold_res))
    best_row <- res[which.max(res$CIndex_1se), ]
    cat(sprintf("done  best alpha=%.1f  CIndex=%.4f\n",
                best_row$alpha, best_row$CIndex_1se)); flush.console()
    res
  }, error=function(e) { cat("error\n"); NULL })
})
en_cx_outer_df  <- do.call(rbind, Filter(Negate(is.null), en_cx_outer))
en_cx_alpha_sum <- aggregate(CIndex_1se ~ alpha, data=en_cx_outer_df, FUN=mean)
best_alpha_cx   <- en_cx_alpha_sum$alpha[which.max(en_cx_alpha_sum$CIndex_1se)]
progress(sprintf("  Cox: best alpha=%.1f", best_alpha_cx))

# Per-fold CV with best_alpha_cx
progress(sprintf("  Cox: per-fold CV with best alpha=%.1f...", best_alpha_cx))

cv_en_cx_raw <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  cat(sprintf("    [EN Cox CV] Fold %d/5 ... ", i)); flush.console()
  tr <- analysis(sp); va <- assessment(sp)
  tryCatch({
    enc    <- encode_predictors(tr, va, "days_to_event")
    X_tr_f <- as.matrix(enc$train); X_va_f <- as.matrix(enc$test)
    y_cx_f <- Surv(tr$days_to_event, tr$lapse_flag)
    set.seed(42 + i)
    inner_clients     <- unique(tr$id)
    client_inner_fold <- setNames(
      sample(rep(1:5, length.out = length(inner_clients))),
      as.character(inner_clients)
    )
    inner_foldid <- client_inner_fold[as.character(tr$id)]
    set.seed(42)
    cv_fit <- cv.glmnet(X_tr_f, y_cx_f, family="cox", alpha=best_alpha_cx,
                        standardize=TRUE, foldid=inner_foldid, type.measure="C")
    make_cox_fold_row <- function(s_val) {
      sel <- intersect(extract_nonzero_vars(cv_fit, s_val),
                       intersect(names(enc$train), names(enc$test)))
      if (length(sel)==0) return(NULL)
      tr_sel <- cbind(time=tr$days_to_event, event=tr$lapse_flag,
                      enc$train[,sel,drop=FALSE])
      va_sel <- cbind(time=va$days_to_event, event=va$lapse_flag,
                      enc$test[,sel,drop=FALSE])
      fold_cox_metrics(tr, va, tr_sel, va_sel)
    }
    res <- list(min=make_cox_fold_row("lambda.min"),
                se1=make_cox_fold_row("lambda.1se"))
    cat(sprintf("done\n")); flush.console()
    res
  }, error=function(e) { cat("error\n"); list(min=NULL, se1=NULL) })
})
cv_raw_en_cx_min <- Filter(Negate(is.null), lapply(cv_en_cx_raw, `[[`, "min"))
cv_raw_en_cx_1se <- Filter(Negate(is.null), lapply(cv_en_cx_raw, `[[`, "se1"))
cv_sum_en_cx_min <- cv_summary(do.call(rbind, cv_raw_en_cx_min))
cv_sum_en_cx_1se <- cv_summary(do.call(rbind, cv_raw_en_cx_1se))

# Full-training Cox with best_alpha_cx
X_cx <- as.matrix(enc_cx$train)
y_cx <- Surv(train_data$days_to_event, train_data$lapse_flag)
set.seed(42)
best_en_cx <- cv.glmnet(X_cx, y_cx, family="cox", alpha=best_alpha_cx,
                         foldid=foldid_cx, type.measure="C", standardize=TRUE)
vars_en_cx_min <- extract_nonzero_vars(best_en_cx, "lambda.min")
vars_en_cx_1se <- extract_nonzero_vars(best_en_cx, "lambda.1se")
orig_en_cx_min <- unique(strip_suffix(vars_en_cx_min, known_cx))
orig_en_cx_1se <- unique(strip_suffix(vars_en_cx_1se, known_cx))
en_cx_min_fit  <- unpen_cox(as.data.frame(train_mod_srv), vars_en_cx_min)
en_cx_1se_fit  <- unpen_cox(as.data.frame(train_mod_srv), vars_en_cx_1se)

make_cox_test_metrics <- function(fit, vars_used) {
  lp_te  <- predict(fit, test_mod_srv, type="lp")
  sp_te  <- get_surv_from_cox(fit, lp_te)
  ci_te  <- cindex_score(test_mod_srv$time, test_mod_srv$event, lp_te)
  ibs_te <- integrated_brier_score(train_mod_srv$time, train_mod_srv$event,
                                   test_mod_srv$time, test_mod_srv$event, fit, lp_te)
  wauc   <- auc_score(test_data$lapse_flag, 1-sp_te$S_T2)
  cond   <- pmin((1-sp_te$S_T1[lapsed_te_idx]) /
                   pmax(1-sp_te$S_T2[lapsed_te_idx],1e-10), 1)
  wnauc  <- auc_score(early_te_obs, cond)
  list(ci=ci_te, ibs=ibs_te, wauc=wauc, wnauc=wnauc,
       n_orig=length(unique(strip_suffix(vars_used, known_cx))))
}
cx_min_te <- make_cox_test_metrics(en_cx_min_fit, vars_en_cx_min)
cx_1se_te <- make_cox_test_metrics(en_cx_1se_fit, vars_en_cx_1se)
row_en_cx_min <- make_row_cox(sprintf("EN a=%.1f lambda.min", best_alpha_cx),
                              cx_min_te$n_orig, cv_sum_en_cx_min,
                              cx_min_te$ci, cx_min_te$ibs, cx_min_te$wauc, cx_min_te$wnauc)
row_en_cx_1se <- make_row_cox(sprintf("EN a=%.1f lambda.1se", best_alpha_cx),
                              cx_1se_te$n_orig, cv_sum_en_cx_1se,
                              cx_1se_te$ci, cx_1se_te$ibs, cx_1se_te$wauc, cx_1se_te$wnauc)
save_results(data.frame(variable=orig_en_cx_min), "elastic_min", "cox", "selected_vars")
save_results(do.call(rbind,cv_raw_en_cx_min),     "elastic_min", "cox", "cv_results")
save_results(row_en_cx_min,                        "elastic_min", "cox", "metrics")
save_results(data.frame(variable=orig_en_cx_1se), "elastic_1se", "cox", "selected_vars")
save_results(do.call(rbind,cv_raw_en_cx_1se),     "elastic_1se", "cox", "cv_results")
save_results(row_en_cx_1se,                        "elastic_1se", "cox", "metrics")
progress(sprintf("  Cox EN complete: lambda.min=%d vars | lambda.1se=%d vars.",
                 length(orig_en_cx_min), length(orig_en_cx_1se)))

# ==============================================================================
# 10.3 Stage 2 Elastic Net — lapsed-only subset
# ==============================================================================
progress("  Stage 2: outer fold x alpha grid (AUC criterion)...")

en_s2_outer <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  cat(sprintf("    [EN S2 alpha search] Fold %d/5 ... ", i)); flush.console()
  tr <- analysis(sp) %>% filter(lapse_flag==1) %>%
        mutate(early_lapse=as.integer(days_to_event<=LANDMARK_T1))
  va <- assessment(sp) %>% filter(lapse_flag==1) %>%
        mutate(early_lapse=as.integer(days_to_event<=LANDMARK_T1))
  if (nrow(tr)<10 || nrow(va)<2 || length(unique(tr$early_lapse))<2) {
    cat("skipped\n"); return(NULL)
  }
  tryCatch({
    enc    <- encode_predictors(tr, va, "early_lapse")
    X_tr_f <- as.matrix(enc$train); X_va_f <- as.matrix(enc$test)
    set.seed(42 + i)
    inner_clients     <- unique(tr$id)
    client_inner_fold <- setNames(
      sample(rep(1:5, length.out = length(inner_clients))),
      as.character(inner_clients)
    )
    inner_foldid <- client_inner_fold[as.character(tr$id)]
    
    fold_res <- lapply(ALPHA_GRID, function(a) {
      set.seed(42)
      cv_fit <- tryCatch(
        cv.glmnet(X_tr_f, tr$early_lapse, family="binomial", alpha=a,
                  standardize=TRUE, foldid=inner_foldid, type.measure="auc"),
        error=function(e) NULL
      )
      if (is.null(cv_fit)) return(NULL)
      p_1se <- as.vector(predict(cv_fit, newx=X_va_f, s="lambda.1se", type="response"))
      data.frame(alpha=a, AUC_1se=auc_score(va$early_lapse, p_1se))
    })
    res      <- do.call(rbind, Filter(Negate(is.null), fold_res))
    best_row <- res[which.max(res$AUC_1se), ]
    cat(sprintf("done  best alpha=%.1f  AUC=%.4f\n",
                best_row$alpha, best_row$AUC_1se)); flush.console()
    res
  }, error=function(e) { cat("error\n"); NULL })
})
en_s2_outer_df  <- do.call(rbind, Filter(Negate(is.null), en_s2_outer))
if (is.null(en_s2_outer_df) || nrow(en_s2_outer_df) == 0) {
  progress("  WARNING: All Stage 2 EN outer folds skipped — defaulting best_alpha_s2 to 0.5.")
  best_alpha_s2 <- 0.5
} else {
  en_s2_alpha_sum <- aggregate(AUC_1se ~ alpha, data=en_s2_outer_df, FUN=mean)
  best_alpha_s2   <- en_s2_alpha_sum$alpha[which.max(en_s2_alpha_sum$AUC_1se)]
}

progress(sprintf("  Stage 2: best alpha=%.1f", best_alpha_s2))

# Per-fold CV with best_alpha_s2
progress(sprintf("  Stage 2: per-fold CV with best alpha=%.1f...", best_alpha_s2))

cv_en_s2_raw <- lapply(seq_along(cv_folds$splits), function(i) {
  sp <- cv_folds$splits[[i]]
  cat(sprintf("    [EN S2 CV] Fold %d/5 ... ", i)); flush.console()
  tr <- analysis(sp) %>% filter(lapse_flag==1) %>%
    mutate(early_lapse=as.integer(days_to_event<=LANDMARK_T1))
  va <- assessment(sp) %>% filter(lapse_flag==1) %>%
    mutate(early_lapse=as.integer(days_to_event<=LANDMARK_T1))
  if (nrow(tr)<10 || nrow(va)<2 || length(unique(tr$early_lapse))<2) {
    cat("skipped\n"); return(list(min=NULL, se1=NULL))
  }
  tryCatch({
    enc    <- encode_predictors(tr, va, "early_lapse")
    X_tr_f <- as.matrix(enc$train); X_va_f <- as.matrix(enc$test)
    set.seed(42 + i)
    inner_clients     <- unique(tr$id)
    client_inner_fold <- setNames(
      sample(rep(1:5, length.out = length(inner_clients))),
      as.character(inner_clients)
    )
    inner_foldid <- client_inner_fold[as.character(tr$id)]
    set.seed(42)
    cv_fit <- cv.glmnet(X_tr_f, tr$early_lapse, family="binomial",
                        alpha=best_alpha_s2, standardize=TRUE, foldid=inner_foldid, type.measure="auc")
    make_fold_row <- function(s_val) {
      sel <- intersect(extract_nonzero_vars(cv_fit, s_val),
                       intersect(names(enc$train), names(enc$test)))
      if (length(sel)==0) return(NULL)
      m  <- glm(y~., data=cbind(y=tr$early_lapse, enc$train[,sel,drop=FALSE]),
                family=binomial())
      p  <- predict(m, enc$test[,sel,drop=FALSE], type="response")
      fm <- binary_lapse_metrics(va$early_lapse, p); fm$AIC <- AIC(m); fm
    }
    res <- list(min=make_fold_row("lambda.min"), se1=make_fold_row("lambda.1se"))
    cat("done\n"); flush.console()
    res
  }, error=function(e) { cat("error\n"); list(min=NULL, se1=NULL) })
})
cv_raw_en_s2_min <- Filter(Negate(is.null), lapply(cv_en_s2_raw, `[[`, "min"))
cv_raw_en_s2_1se <- Filter(Negate(is.null), lapply(cv_en_s2_raw, `[[`, "se1"))
cv_sum_en_s2_min <- if (length(cv_raw_en_s2_min)>0) cv_summary(do.call(rbind,cv_raw_en_s2_min)) else NULL
cv_sum_en_s2_1se <- if (length(cv_raw_en_s2_1se)>0) cv_summary(do.call(rbind,cv_raw_en_s2_1se)) else NULL

# Full-training Stage 2 with best_alpha_s2
X_s2 <- as.matrix(enc_s2$train); y_s2 <- train_data_lapsed$early_lapse
set.seed(42)
best_en_s2 <- cv.glmnet(X_s2, y_s2, family="binomial", alpha=best_alpha_s2,
                         foldid=foldid_s2, type.measure="auc", standardize=TRUE)
vars_en_s2_min <- extract_nonzero_vars(best_en_s2, "lambda.min")
vars_en_s2_1se <- extract_nonzero_vars(best_en_s2, "lambda.1se")
orig_en_s2_min <- unique(strip_suffix(vars_en_s2_min, known_s2))
orig_en_s2_1se <- unique(strip_suffix(vars_en_s2_1se, known_s2))

make_en_s2_row <- function(label, vars, orig_vars, cv_sum_obj) {
  if (length(vars)==0 || is.null(cv_sum_obj)) {
    return(data.frame(Method=label, N_Vars=0,
                      CV_AUC_Mean=NA, CV_AUC_SE=NA,
                      CV_Brier_Mean=NA, CV_Brier_SE=NA,
                      CV_F1_Mean=NA, CV_F1_SE=NA,
                      Test_AUC=NA, Test_Brier=NA, Test_F1=NA,
                      stringsAsFactors=FALSE))
  }
  fit      <- unpen_glm(as.data.frame(train_s2), vars, "y")
  test_res <- test_eval_binary(fit, test_s2, "y")
  make_row_binary(label, length(orig_vars), cv_sum_obj, test_res)
}

row_en_s2_min <- make_en_s2_row(sprintf("EN a=%.1f lambda.min", best_alpha_s2),
                                vars_en_s2_min, orig_en_s2_min, cv_sum_en_s2_min)
row_en_s2_1se <- make_en_s2_row(sprintf("EN a=%.1f lambda.1se", best_alpha_s2),
                                vars_en_s2_1se, orig_en_s2_1se, cv_sum_en_s2_1se)

save_results(data.frame(variable=orig_en_s2_min), "elastic_min", "stage2", "selected_vars")
if (length(cv_raw_en_s2_min)>0)
  save_results(do.call(rbind,cv_raw_en_s2_min), "elastic_min", "stage2", "cv_results")
save_results(row_en_s2_min,                      "elastic_min", "stage2", "metrics")
save_results(data.frame(variable=orig_en_s2_1se), "elastic_1se", "stage2", "selected_vars")
if (length(cv_raw_en_s2_1se)>0)
  save_results(do.call(rbind,cv_raw_en_s2_1se), "elastic_1se", "stage2", "cv_results")
save_results(row_en_s2_1se,                      "elastic_1se", "stage2", "metrics")
progress(sprintf("  Stage 2 EN complete: lambda.min=%d vars | lambda.1se=%d vars.",
                 length(orig_en_s2_min), length(orig_en_s2_1se)))

# ==============================================================================
# SECTION 11: COMPARISON TABLES
# Best value in each column is highlighted green. Higher-is-better: AUC, C-Index,
# W-AUC, WN-AUC, F1. Lower-is-better: Brier, IBS, and all SE columns. N_Vars is
# left unhighlighted (parsimony is handled separately by the 1-SE rule below).
# ==============================================================================
progress("=== BUILDING COMPARISON TABLES ===")

# ------------------------------------------------------------------------------
# Helper: apply green highlighting to the best cell in each nominated column
# ------------------------------------------------------------------------------
highlight_best <- function(df, higher_cols, lower_cols, se_cols) {
  
  apply_hl <- function(df_out, col, best_fn) {
    if (!col %in% names(df_out)) return(df_out)
    raw  <- suppressWarnings(as.numeric(df_out[[col]]))
    idx  <- best_fn(raw)                     # which.max / which.min, NA-aware
    if (length(idx) == 0 || all(is.na(raw))) return(df_out)
    fmt  <- ifelse(is.na(raw), "—", formatC(raw, digits = 4, format = "f"))
    df_out[[col]] <- cell_spec(
      fmt,
      format     = "html",
      background = ifelse(seq_len(nrow(df_out)) == idx, "#27ae60", "transparent"),
      color      = ifelse(seq_len(nrow(df_out)) == idx, "white",   "inherit"),
      bold       = seq_len(nrow(df_out)) == idx
    )
    df_out
  }
  
  df_hl <- df
  for (col in higher_cols)
    df_hl <- apply_hl(df_hl, col, function(v) which.max(replace(v, is.na(v), -Inf)))
  for (col in c(lower_cols, se_cols))
    df_hl <- apply_hl(df_hl, col, function(v) which.min(replace(v, is.na(v),  Inf)))
  df_hl
}

# ------------------------------------------------------------------------------
# Updated save_comparison_table — accepts metric direction vectors
# ------------------------------------------------------------------------------
save_comparison_table <- function(comparison_df, model_name, col_labels,
                                  higher_cols, lower_cols, se_cols) {
  csv_path  <- here(sprintf("comparison_%s.csv",  model_name))
  rds_path  <- here(sprintf("comparison_%s.rds",  model_name))
  html_path <- here(sprintf("comparison_%s.html", model_name))
  
  # Save raw (un-highlighted) data
  write.csv(comparison_df, csv_path, row.names = FALSE)
  saveRDS(comparison_df,   rds_path)
  
  # Build highlighted display copy
  df_hl <- highlight_best(comparison_df, higher_cols, lower_cols, se_cols)
  
  tbl_html <- kable(
    df_hl,
    format     = "html",
    escape     = FALSE,          # <-- needed so cell_spec HTML renders
    col.names  = col_labels,
    caption    = sprintf("Feature Selection Comparison — %s", model_name),
    digits     = 4,
    align      = c("l", rep("c", ncol(df_hl) - 1))
  ) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE, font_size = 13
    ) %>%
    column_spec(1, bold = TRUE) %>%
    row_spec(0, bold = TRUE, background = "#2c3e50", color = "white")
  
  save_kable(tbl_html, html_path)
  progress(sprintf("  Saved: comparison_%s.csv / .rds / .html", model_name))
  invisible(comparison_df)
}

# ==============================================================================
# Stage 1 — P(lapse = 1)
# ==============================================================================
s1_higher <- c("CV_AUC_Mean",   "CV_F1_Mean",   "Test_AUC",   "Test_F1")
s1_lower  <- c("CV_Brier_Mean", "Test_Brier")
s1_se     <- c("CV_AUC_SE",     "CV_Brier_SE",  "CV_F1_SE")

comp_s1 <- rbind(row_base_s1, row_step_s1, row_bor_s1, row_en_s1_min, row_en_s1_1se)
comp_s1_labels <- c(
  "Method", "N Vars",
  "CV AUC Mean",   "CV AUC SE",
  "CV Brier Mean", "CV Brier SE",
  "CV F1 Mean",    "CV F1 SE",
  "Test AUC", "Test Brier", "Test F1"
)
save_comparison_table(comp_s1, "stage1", comp_s1_labels, s1_higher, s1_lower, s1_se)

cat("\n--- Stage 1 Feature Selection Comparison ---\n")
print(kable(comp_s1, col.names = comp_s1_labels, digits = 4,
            caption = "Stage 1: P(lapse=1) — Feature Selection Comparison"))

# ==============================================================================
# Cox PH — days_to_event
# ==============================================================================
cx_higher <- c("CV_CIndex_Mean", "CV_W_AUC_Mean", "CV_WN_AUC_Mean",
               "Test_CIndex",    "Test_W_AUC",    "Test_WN_AUC")
cx_lower  <- c("CV_IBS_Mean",    "Test_IBS")
cx_se     <- c("CV_CIndex_SE",   "CV_IBS_SE",     "CV_W_AUC_SE",  "CV_WN_AUC_SE")

comp_cx <- rbind(row_base_cx, row_step_cx, row_bor_cx, row_en_cx_min, row_en_cx_1se)
comp_cx_labels <- c(
  "Method", "N Vars",
  "CV C-Index Mean", "CV C-Index SE",
  "CV IBS Mean",     "CV IBS SE",
  "CV W-AUC Mean",   "CV W-AUC SE",
  "CV WN-AUC Mean",  "CV WN-AUC SE",
  "Test C-Index", "Test IBS", "Test W-AUC", "Test WN-AUC"
)
save_comparison_table(comp_cx, "cox", comp_cx_labels, cx_higher, cx_lower, cx_se)

cat("\n--- Cox PH Feature Selection Comparison ---\n")
print(kable(comp_cx, col.names = comp_cx_labels, digits = 4,
            caption = "Cox PH — Feature Selection Comparison"))

# ==============================================================================
# Stage 2 — P(early lapse | lapse = 1)
# ==============================================================================
s2_higher <- c("CV_AUC_Mean",   "CV_F1_Mean",   "Test_AUC",   "Test_F1")
s2_lower  <- c("CV_Brier_Mean", "Test_Brier")
s2_se     <- c("CV_AUC_SE",     "CV_Brier_SE",  "CV_F1_SE")

comp_s2 <- rbind(row_base_s2, row_step_s2, row_bor_s2, row_en_s2_min, row_en_s2_1se)
comp_s2_labels <- c(
  "Method", "N Vars",
  "CV AUC Mean",   "CV AUC SE",
  "CV Brier Mean", "CV Brier SE",
  "CV F1 Mean",    "CV F1 SE",
  "Test AUC", "Test Brier", "Test F1"
)
save_comparison_table(comp_s2, "stage2", comp_s2_labels, s2_higher, s2_lower, s2_se)

cat("\n--- Stage 2 Feature Selection Comparison ---\n")
print(kable(comp_s2, col.names = comp_s2_labels, digits = 4,
            caption = "Stage 2: P(early lapse | lapse=1) — Feature Selection Comparison"))

# ==============================================================================
# SECTION 12: PREFERRED MODEL SELECTION (1-SE Rule)
# ==============================================================================

progress("=== PREFERRED MODEL SELECTION (1-SE Rule) ===")

select_preferred <- function(comp_df, metric_mean_col, metric_se_col,
                             higher_is_better = TRUE) {
  vals <- comp_df[[metric_mean_col]]
  ses  <- comp_df[[metric_se_col]]
  if (higher_is_better) {
    best_val  <- max(vals, na.rm=TRUE)
    threshold <- best_val - ses[which.max(vals)]
    eligible  <- which(vals >= threshold)
  } else {
    best_val  <- min(vals, na.rm=TRUE)
    threshold <- best_val + ses[which.min(vals)]
    eligible  <- which(vals <= threshold)
  }
  preferred <- eligible[which.min(comp_df$N_Vars[eligible])]
  comp_df$Method[preferred]
}

pref_s1 <- select_preferred(comp_s1, "CV_AUC_Mean",    "CV_AUC_SE")
pref_cx <- select_preferred(comp_cx, "CV_CIndex_Mean",  "CV_CIndex_SE")
pref_s2 <- select_preferred(comp_s2, "CV_AUC_Mean",    "CV_AUC_SE")

preferred_summary <- data.frame(
  Model     = c("Stage 1 — P(lapse=1)",
                "Cox PH — days_to_event",
                "Stage 2 — P(early lapse | lapse=1)"),
  Preferred = c(pref_s1, pref_cx, pref_s2),
  stringsAsFactors = FALSE
)

write.csv(preferred_summary, here("fs_preferred_models.csv"), row.names=FALSE)
saveRDS(preferred_summary,  here("fs_preferred_models.rds"))
cat("\n--- Preferred Model per Component (1-SE Rule) ---\n")
print(preferred_summary)

progress("Feature selection script complete. All results saved.")
