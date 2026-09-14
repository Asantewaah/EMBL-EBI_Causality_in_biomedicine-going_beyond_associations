# utils.R
#
# Single source of truth for:
#   • package loading
#   • scalar math helpers      (.logit, .expit, .scale_y, .pvalue, .clip_ps)
#   • data helpers             (.W_matrix, .make_QX)
#   • shared NA result builder (.na_est)
#   • collaborative CV loop    (.collab_cv_loss)
#   • H-tilde covariate        (.make_H_tilde)
#   • cross-fitting            (.crossfit_Q)
#   • cross-validated target   (.cv_target)
#   • H-tilde targeting loop   (.cv_target_htilde)
#   • IC-based inference       (.ic_inference)
#   • true ATE calculator      (compute_true_ate)
#   • result collection        (collect_results_snp, collect_results)
#   • MC summary               (MC_res_anchorgene)
#   • plotting helpers
#
# Every estimator file does: source("utils.R")
# DO NOT duplicate package calls or any of these helpers elsewhere.


suppressPackageStartupMessages({
  library(glmnet)
  library(glm2)
  library(caret)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ctmle)
})


# Scalar math helpers


.logit <- function(p) {
  p <- pmax(pmin(p, 0.999), 0.001)
  log(p / (1 - p))
}

.expit <- function(x) 1 / (1 + exp(-x))

.scale_y <- function(y) {
  a <- min(y); b <- max(y)
  list(y_scaled = (y - a) / (b - a), a = a, b = b)
}

.pvalue <- function(est, se) {
  if (is.na(se) || se <= 0) return(NA_real_)
  2 * pnorm(-abs(est / se))
}

# Clip propensity scores to [lo, hi] (default [0.01, 0.99])
.clip_ps <- function(x, lo = 0.01, hi = 0.99) pmax(pmin(x, hi), lo)


# Data helpers


# Extract the W covariate matrix (columns whose names start with "W")
.W_matrix <- function(data) {
  w_cols <- grep("^W", colnames(data), value = TRUE)
  as.matrix(data[, w_cols, drop = FALSE])
}

# Build the three model matrices needed by every LASSO-outcome estimator.
# Returns a list(QX, QX0, QX1).
.make_QX <- function(A, X) {
  list(
    QX  = model.matrix(~ A + ., data = data.frame(A = A, X)),
    QX0 = model.matrix(~ A + ., data = data.frame(A = 0, X)),
    QX1 = model.matrix(~ A + ., data = data.frame(A = 1, X))
  )
}


# NA-result builder
#
# Usage: .na_est("C-TMLE1", best_lambda = NA_real_, n_forced = NA_integer_)
# Extra named args are appended to the list as-is.

.na_est <- function(estimator, ...) {
  base <- list(
    est      = NA_real_,
    se       = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    pvalue   = NA_real_,
    estimator = estimator
  )
  extras <- list(...)
  c(base, extras)
}

# .assign_ld_blocks()
#
# Partitions SNPs into LD blocks using single-linkage on adjacent r².
# SNP j joins the current block if r²(j-1, j) >= threshold.
#
# Arguments:
#   r2_mat    : p x p matrix, rows/cols in genomic order
#   snp_ids   : character vector of SNP IDs (length p), used as names
#   threshold : r² cutoff (default 0.1)
#
# Returns: named integer vector of length p

# .assign_ld_blocks()
.assign_ld_blocks <- function(r2_mat, snp_ids, threshold = 0.1) {
  p           <- nrow(r2_mat)
  block_id    <- integer(p)
  current     <- 1L
  block_id[1] <- current

  # Guard: 2:1 in R counts DOWN — must check p > 1 explicitly
  if (p > 1) {
    for (j in 2:p) {
      if (r2_mat[j - 1, j] >= threshold) {
        block_id[j] <- current
      } else {
        current     <- current + 1L
        block_id[j] <- current
      }
    }
  }

  names(block_id) <- snp_ids
  block_id
}


# .load_ld_blocks()
#
# Loads r2_mat.csv and snp_info.csv, subsets to a SNP window,
# and returns a named integer vector of LD block IDs.
#
# Matching is positional — both files are in the same SNP order
# as written by simulate_genotypes_CEU.R.
#
# Arguments:
#   r2_path       : path to r2_mat.csv
#   snp_info_path : path to snp_info.csv
#   snp_window    : integer vector of SNP indices to subset (e.g. 201:400)
#                   NULL means use all SNPs
#   threshold     : r² cutoff for block merging (default 0.1)

.load_ld_blocks <- function(r2_path, snp_info_path,
                            snp_window = NULL, threshold = 0.1) {

  snp_info <- read.csv(snp_info_path, row.names = 1,
                       stringsAsFactors = FALSE)

  # r2_mat.csv has duplicate kb-position row labels when multiple SNPs
  # round to the same kb value.  read.csv(..., row.names = 1) errors on
  # duplicates, so we read without row names and drop the label column.
  r2_raw <- read.csv(r2_path, header = TRUE,
                     row.names = NULL,      # <-- key change
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

  # First column is the kb label — drop it, keep the numeric matrix
  r2_mat <- as.matrix(r2_raw[, -1, drop = FALSE])

  if (nrow(r2_mat) != nrow(snp_info))
    stop(".load_ld_blocks: r2_mat has ", nrow(r2_mat), " rows but ",
         "snp_info has ", nrow(snp_info), " rows.")

  if (!is.null(snp_window)) {
    r2_mat   <- r2_mat[snp_window, snp_window, drop = FALSE]
    snp_info <- snp_info[snp_window, , drop = FALSE]
  }

  .assign_ld_blocks(r2_mat,
                    snp_ids   = snp_info$snp_id,
                    threshold = threshold)
}



# .ic_inference()
#
# Efficient influence function (EIF) -based SE and CI.
# mean_IF should be ~0 after a successful targeting step; large values
# indicate non-convergence or a numerical issue.

.ic_inference <- function(A, Y, Q_star1, Q_star0, g1W,
                          ATE_n, b_minus_a, alpha) {      

  IF  <- (A / g1W) * (Y - Q_star1) -
         ((1 - A) / (1 - g1W)) * (Y - Q_star0) +
         (Q_star1 - Q_star0) - ATE_n

  SE <- sqrt(var(IF) / length(A)) * b_minus_a

  z   <- qnorm(1 - alpha / 2)
  ATE <- ATE_n * b_minus_a
  list(
    est      = ATE,
    se       = SE,
    ci_lower = ATE - z * SE,
    ci_upper = ATE + z * SE,
    pvalue   = .pvalue(ATE, SE),
    mean_IF  = mean(IF)
  )
}

.ic_inference_cv <- function(A, Y, Q_star1, Q_star0, g1W, folds,
                              ATE_n, b_minus_a, alpha) {

  n  <- length(A)
  IF <- numeric(n)

  for (v in seq_along(folds)) {
    val <- folds[[v]]
    IF[val] <- (A[val] / g1W[val]) * (Y[val] - Q_star1[val]) -
               ((1 - A[val]) / (1 - g1W[val])) * (Y[val] - Q_star0[val]) +
               (Q_star1[val] - Q_star0[val]) - ATE_n
  }

  SE <- sqrt(mean(IF^2) / n) * b_minus_a   # CV_Var formula

  z   <- qnorm(1 - alpha / 2)
  ATE <- ATE_n * b_minus_a
  list(
    est      = ATE,
    se       = SE,
    ci_lower = ATE - z * SE,
    ci_upper = ATE + z * SE,
    pvalue   = .pvalue(ATE, SE),
    mean_IF  = mean(IF)
  )
}


# .collab_cv_loss()
#
# Inner loop shared by ctmle1(), ctmle0star(), and ctmle1_oal().
#
# For each fold: fits a propensity model at a single lambda on the training
# set, runs a one-step targeting update on the validation set, and returns
# the squared-error loss.  The caller iterates over lambda_seq and picks the
# lambda that minimises the accumulated loss.
#
# Arguments:
#   X          : covariate matrix [n x p]
#   A, Y       : treatment and scaled outcome [length n]
#   Q1W, Q0W   : cross-fitted Q predictions [length n]
#   folds      : fold list (list of held-out index vectors)
#   lambda_seq : numeric vector of lambda values to evaluate
#   ps_fit_fn  : function(X_tr, A_tr, lambda) -> glmnet model
#                Default fits a standard LASSO logistic regression.
#                ctmle1_oal() passes a closure that includes penalty.factor.
#
# Returns: best_lambda (scalar)

.collab_cv_loss <- function(X, A, Y, Q1W, Q0W, folds,
                            lambda_seq,
                            ps_fit_fn = function(X_tr, A_tr, lam)
                              glmnet(X_tr, A_tr, family = "binomial",
                                     lambda = lam)) {
  n          <- nrow(X)
  V          <- length(folds)
  best_loss  <- Inf
  best_lambda <- lambda_seq[length(lambda_seq)]   # fallback = smallest lambda

  for (lam in lambda_seq) {
    cv_loss     <- 0
    valid_folds <- 0L

    for (v in seq_len(V)) {
      tr  <- setdiff(seq_len(n), folds[[v]])
      val <- folds[[v]]

      g_fit <- tryCatch(ps_fit_fn(X[tr, , drop = FALSE], A[tr], lam),
                        error = function(e) NULL)
      if (is.null(g_fit)) next

      g1W_val <- .clip_ps(
        as.numeric(predict(g_fit, X[val, , drop = FALSE], type = "response"))
      )

      Q1_val <- Q1W[val]; Q0_val <- Q0W[val]

      d1 <- data.frame(Yv = Y[val], H1 = A[val],       off1 = .logit(Q1_val))
      d0 <- data.frame(Yv = Y[val], H0 = 1 - A[val],   off0 = .logit(Q0_val))

      fit1 <- tryCatch(
        glm2(Yv ~ -1 + H1 + offset(off1), family = quasibinomial,
             data = d1, weights = 1 / g1W_val),
        error = function(e) NULL)
      fit0 <- tryCatch(
        glm2(Yv ~ -1 + H0 + offset(off0), family = quasibinomial,
             data = d0, weights = 1 / (1 - g1W_val)),
        error = function(e) NULL)

      if (is.null(fit1) || is.null(fit0)) next
      eps1 <- coef(fit1)[1]; eps0 <- coef(fit0)[1]
      if (!is.finite(eps1) || !is.finite(eps0)) next

      Qs1 <- .expit(.logit(Q1_val) + eps1)
      Qs0 <- .expit(.logit(Q0_val) + eps0)
      fold_loss <- mean((Y[val] - (A[val] * Qs1 + (1 - A[val]) * Qs0))^2)
      if (!is.finite(fold_loss)) next

      cv_loss     <- cv_loss + fold_loss
      valid_folds <- valid_folds + 1L
    }

    if (valid_folds == 0L) next
    if (cv_loss < best_loss) {
      best_loss   <- cv_loss
      best_lambda <- lam
    }
  }

  best_lambda
}


# .make_H_tilde()
#
# The "extra clever covariate" used by ctmle0() and ctmle0star().
#
# Arguments:
#   A      : treatment vector
#   g1W    : propensity at best_lambda
#   ps_fit : a cv.glmnet propensity fit (used to predict at lam_d)
#   lam    : the selected lambda (log-scale nudge computed internally)
#
# Returns: numeric vector of length n

.make_H_tilde <- function(A, g1W, ps_fit, lam) {
  lam_d <- exp(log(lam) + 1e-3)
  lam_d <- min(lam_d, max(ps_fit$glmnet.fit$lambda))
  g1W_d <- .clip_ps(
    as.numeric(predict(ps_fit, newx = attr(g1W, "X"), s = lam_d,
                       type = "response"))
  )
  dg <- g1W_d - g1W
  ((1 - A) / (1 - g1W)^2) * (-dg) + (A / g1W^2) * dg
}

# Note: .make_H_tilde() requires X to be attached as an attribute of g1W by
# the caller (see ctmle0 / ctmle0star).  Alternatively callers can compute
# H_tilde directly using the two-liner below — both approaches are equivalent:
#
#   lam_d   <- min(exp(log(lam) + 1e-3), max(ps_fit$glmnet.fit$lambda))
#   g1W_d   <- .clip_ps(predict(ps_fit, newx = X, s = lam_d, type = "response"))
#   H_tilde <- ((1 - A) / (1 - g1W)^2) * (g1W_d - g1W) * c(-1, 1)[A + 1]
#              # simplified form; the explicit version is in each estimator.


# .crossfit_Q()
#
# Cross-fitted outcome model predictions.
#
# Fits cv.glmnet on (V-1) training folds and predicts Q(A=0,W) and Q(A=1,W)
# on the held-out fold, ensuring the EIF residuals Y - Q_star are
# out-of-sample.  Also returns Q_full (in-sample, full-data fit) used only
# for residual ordering in scalable_ctmle().
#
# Arguments:
#   QX, QX0, QX1 : model matrices [n x k] (from .make_QX())
#   Y            : scaled outcome [length n]
#   V            : cross-fitting folds (default 10)
#   folds        : optional pre-built fold list (reused when provided)
#
# Returns: list(Q0W, Q1W, Q_full, folds)

.crossfit_Q <- function(QX, QX0, QX1, Y, V = 10, folds = NULL) {
  n   <- nrow(QX)
  Q0W <- numeric(n)
  Q1W <- numeric(n)

  if (is.null(folds))
    folds <- caret::createFolds(Y, k = V, list = TRUE, returnTrain = FALSE)

  for (v in seq_len(V)) {
    tr  <- setdiff(seq_len(n), folds[[v]])
    val <- folds[[v]]
    Qf  <- tryCatch(
      cv.glmnet(QX[tr, , drop = FALSE], Y[tr], alpha = 1, nfolds = 5),
      error = function(e) NULL
    )
    if (is.null(Qf)) next
    Q0W[val] <- as.numeric(predict(Qf, newx = QX0[val, , drop = FALSE], s = "lambda.min"))
    Q1W[val] <- as.numeric(predict(Qf, newx = QX1[val, , drop = FALSE], s = "lambda.min"))
  }

  # Full-data fit: used only for residual-ordering in scalable_ctmle()
  Q_init <- cv.glmnet(QX, Y, alpha = 1, nfolds = V)
  Q_full <- as.numeric(predict(Q_init, newx = QX, s = "lambda.min"))

  list(Q0W = Q0W, Q1W = Q1W, Q_full = Q_full, folds = folds)
}


# .auc()
#
# Rank-based (Mann-Whitney) AUC for a binary label vs. a continuous score.
# Used to bound the OAL undersmoothing search in ctmle1_oal() (Wyss et al.
# 2024, Appendix S1): candidate lambdas are only considered while the
# same-sample AUC of the treatment model stays <= 0.999, the cutoff the
# paper uses as a proxy for stochastic positivity violation.
#
# Arguments:
#   y : binary label (0/1)
#   p : predicted score (e.g. predicted propensity)
#
# Returns: scalar AUC, or NA if either class is empty

.auc <- function(y, p) {
  pos <- p[y == 1]; neg <- p[y == 0]
  n1  <- length(pos); n0 <- length(neg)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}


# .crossfit_g()
#
# Out-of-fold ("cross-fitted") propensity score predictions at a single
# lambda, using a caller-supplied glmnet-fitting closure. Mirrors
# .crossfit_Q() but for the treatment model: fits on the training folds,
# predicts on the held-out fold, so the propensity score for an observation
# never comes from a model that has already seen that observation.
#
# This is the piece that distinguishes Wyss et al. (2024)'s "CF" (cross-fit)
# models (5-8) from their non-CF counterparts (1-4). Their Table 1 shows
# the distinction is not cosmetic: the non-CF collaborative-controlled OAL
# (Model 4) has 0-1.6% coverage, while the cross-fit version (Model 8) has
# 94-97% coverage — cross-fitting the undersmoothed propensity model is
# what prevents the nonoverlap that undersmoothing would otherwise induce.
#
# Arguments:
#   X, A      : covariate matrix and treatment vector
#   folds     : fold list — reuse the same folds used elsewhere in the
#               caller (e.g. the folds passed to .collab_cv_loss()) so g
#               and any other cross-fit nuisance share a single split
#   lam       : a single lambda value
#   ps_fit_fn : function(X_tr, A_tr, lam) -> glmnet model. Default fits a
#               plain LASSO logistic regression; callers needing a
#               penalty.factor (e.g. ctmle1_oal()) pass a closure.
#
# Returns: numeric vector of length n, clipped to [0.01, 0.99]

.crossfit_g <- function(X, A, folds, lam,
                        ps_fit_fn = function(X_tr, A_tr, lam)
                          glmnet(X_tr, A_tr, family = "binomial", lambda = lam)) {
  n   <- nrow(X)
  g1W <- numeric(n)

  for (v in seq_along(folds)) {
    tr  <- setdiff(seq_len(n), folds[[v]])
    val <- folds[[v]]

    g_fit <- tryCatch(ps_fit_fn(X[tr, , drop = FALSE], A[tr], lam),
                      error = function(e) NULL)
    if (is.null(g_fit)) {
      g1W[val] <- mean(A[tr])
      next
    }
    g1W[val] <- as.numeric(predict(g_fit, X[val, , drop = FALSE], type = "response"))
  }

  .clip_ps(g1W)
}

# .crossfit_g_lm()
#
# Out-of-fold propensity predictions for a plain OLS propensity model on a
# FIXED variable set, i.e. Wang-CTMLE's parameterization (A ~ X[, g_vars]
# via lm(), not a glmnet lambda path). Sibling to .crossfit_g(), which
# assumes predict.glmnet's newx/type="response" interface and doesn't fit
# lm()'s predict(newdata = data.frame(...)) signature.
#
# Wang-CTMLE selects g_vars via a greedy forward search + CV vote BEFORE
# this is called; this function does NOT reselect variables, it only
# removes the in-sample-fit/in-sample-predict overfitting in the WINNING
# model's fitted propensity scores. Selection-variability (which variables
# get chosen) is a separate, unresolved caveat -- see wang_ctmle()'s
# boot_var block.
#
# Arguments:
#   X, A   : covariate matrix and treatment vector
#   folds  : fold list (reuse the SAME folds used elsewhere, e.g. folds_cf)
#   g_vars : character vector of column names in X to include (may be
#            character(0) for the intercept-only model)
#
# Returns: numeric vector of length n. NOT passed through .clip_ps() --
# Wang-CTMLE's gW is not clipped anywhere else in this file either; matches
# existing convention.

.crossfit_g_lm <- function(X, A, folds, g_vars) {
  n  <- nrow(X)
  gW <- numeric(n)

  for (v in seq_along(folds)) {
    tr  <- setdiff(seq_len(n), folds[[v]])
    val <- folds[[v]]

    if (length(g_vars) == 0) {
      gW[val] <- mean(A[tr])
      next
    }

    gf <- tryCatch(lm(A[tr] ~ X[tr, g_vars, drop = FALSE]),
                   error = function(e) NULL)
    if (is.null(gf)) {
      gW[val] <- mean(A[tr])
      next
    }
    pred <- tryCatch(
      as.numeric(cbind(1, X[val, g_vars, drop = FALSE]) %*% coef(gf)),
      error = function(e) NULL)
    gW[val] <- if (is.null(pred) || anyNA(pred)) mean(A[tr]) else pred
  }

  gW
}


# .cv_target()
#
# Cross-validated targeting step (simple: no H_tilde).
#
# For each fold: fits epsilon on the training observations (using the
# cross-fitted Q as offset), then applies it to the held-out observations.
# Q_star is therefore built from entirely out-of-fold epsilons.
#
# Arguments:
#   Y, A         : scaled outcome and treatment
#   Q1W, Q0W    : cross-fitted Q predictions (from .crossfit_Q())
#   g1W          : propensity scores [length n]
#   folds        : fold list (must be the same splits used in .crossfit_Q())
#
# Returns: list(Q_star1, Q_star0)

.cv_target <- function(Y, A, Q1W, Q0W, g1W, folds) {
  n       <- length(Y)
  V       <- length(folds)
  Q_star1 <- numeric(n)
  Q_star0 <- numeric(n)

  for (v in seq_len(V)) {
    tr  <- setdiff(seq_len(n), folds[[v]])
    val <- folds[[v]]

    g1W_tr <- .clip_ps(g1W[tr])

    d1 <- data.frame(Ytr = Y[tr], H1 = A[tr],     off1 = .logit(Q1W[tr]))
    d0 <- data.frame(Ytr = Y[tr], H0 = 1 - A[tr], off0 = .logit(Q0W[tr]))

    fit1 <- tryCatch(
      glm2(Ytr ~ -1 + H1 + offset(off1), family = quasibinomial,
           data = d1, weights = 1 / g1W_tr),
      error = function(e) NULL)
    fit0 <- tryCatch(
      glm2(Ytr ~ -1 + H0 + offset(off0), family = quasibinomial,
           data = d0, weights = 1 / (1 - g1W_tr)),
      error = function(e) NULL)

    if (is.null(fit1) || is.null(fit0)) {
      Q_star1[val] <- Q1W[val]
      Q_star0[val] <- Q0W[val]
      next
    }

    eps1 <- coef(fit1)[1]; if (!is.finite(eps1)) eps1 <- 0
    eps0 <- coef(fit0)[1]; if (!is.finite(eps0)) eps0 <- 0

    Q_star1[val] <- .expit(.logit(Q1W[val]) + eps1)
    Q_star0[val] <- .expit(.logit(Q0W[val]) + eps0)
  }

  list(Q_star1 = Q_star1, Q_star0 = Q_star0)
}


# .cv_target_htilde()
#
# Cross-validated targeting step WITH H_tilde (used by ctmle0 / ctmle0star).
#
# Identical to .cv_target() except the GLM includes an extra covariate H_tilde
# on both treatment and control arms, capturing the gradient of the propensity
# loss with respect to lambda.
#
# Arguments:
#   Y, A         : scaled outcome and treatment
#   Q1W, Q0W    : cross-fitted Q predictions
#   g1W          : propensity scores [length n]
#   H_tilde      : H_tilde vector [length n] (from .make_H_tilde() or inline)
#   folds        : fold list
#
# Returns: list(Q_star1, Q_star0)

.cv_target_htilde <- function(Y, A, Q1W, Q0W, g1W, H_tilde, folds) {
  n       <- length(Y)
  V       <- length(folds)
  Q_star1 <- numeric(n)
  Q_star0 <- numeric(n)

  for (v in seq_len(V)) {
    tr  <- setdiff(seq_len(n), folds[[v]])
    val <- folds[[v]]

    g1W_tr <- .clip_ps(g1W[tr])
    Ht_tr  <- H_tilde[tr]
    Ht_val <- H_tilde[val]

    d1 <- data.frame(Ytr = Y[tr], H1 = A[tr],
                     Ht = Ht_tr, off = .logit(Q1W[tr]))
    d0 <- data.frame(Ytr = Y[tr], H0 = 1 - A[tr],
                     Ht = Ht_tr, off = .logit(Q0W[tr]))

    fit1 <- tryCatch(
      glm2(Ytr ~ -1 + H1 + Ht + offset(off), family = quasibinomial,
           data = d1, weights = 1 / g1W_tr),
      error = function(e) NULL)
    fit0 <- tryCatch(
      glm2(Ytr ~ -1 + H0 + Ht + offset(off), family = quasibinomial,
           data = d0, weights = 1 / (1 - g1W_tr)),
      error = function(e) NULL)

    if (is.null(fit1) || is.null(fit0)) {
      Q_star1[val] <- Q1W[val]
      Q_star0[val] <- Q0W[val]
      next
    }

    eps1  <- coef(fit1)[1]; eps1g <- coef(fit1)[2]
    eps0  <- coef(fit0)[1]; eps0g <- coef(fit0)[2]
    if (!is.finite(eps1))  eps1  <- 0
    if (!is.finite(eps1g)) eps1g <- 0
    if (!is.finite(eps0))  eps0  <- 0
    if (!is.finite(eps0g)) eps0g <- 0

    Q_star1[val] <- .expit(.logit(Q1W[val]) + eps1 + eps1g * Ht_val)
    Q_star0[val] <- .expit(.logit(Q0W[val]) + eps0 + eps0g * Ht_val)
  }

  list(Q_star1 = Q_star1, Q_star0 = Q_star0)
}


# .final_target()
#
# Final (full-data) targeting step shared by ctmle1(), ctmle1_oal(), and
# scalable_ctmle().  Fits one epsilon per arm on the full data and returns
# Q_star1, Q_star0.
#
# Arguments:
#   Y, A     : scaled outcome and treatment
#   Q1W, Q0W : cross-fitted Q predictions
#   g1W      : propensity scores [length n]
#
# Returns: list(Q_star1, Q_star0)

.final_target <- function(Y, A, Q1W, Q0W, g1W) {
  g1W <- .clip_ps(g1W)

  d1f <- data.frame(Y = Y, H1 = A,     off1 = .logit(Q1W))
  d0f <- data.frame(Y = Y, H0 = 1 - A, off0 = .logit(Q0W))

  fit1f <- tryCatch(
    glm2(Y ~ -1 + H1 + offset(off1), family = quasibinomial,
         data = d1f, weights = 1 / g1W),
    error = function(e) NULL)
  fit0f <- tryCatch(
    glm2(Y ~ -1 + H0 + offset(off0), family = quasibinomial,
         data = d0f, weights = 1 / (1 - g1W)),
    error = function(e) NULL)

  if (is.null(fit1f) || is.null(fit0f)) return(NULL)

  eps1f <- coef(fit1f)[1]; if (!is.finite(eps1f)) eps1f <- 0
  eps0f <- coef(fit0f)[1]; if (!is.finite(eps0f)) eps0f <- 0

  list(
    Q_star1 = .expit(.logit(Q1W) + eps1f),
    Q_star0 = .expit(.logit(Q0W) + eps0f)
  )
}


# .final_target_htilde()
#
# Full-data (non-CV) targeting step WITH H_tilde — companion to
# .final_target(), mirroring how .cv_target_htilde() relates to .cv_target().
# Needed because ctmle0()/ctmle0star() target with an extra H_tilde covariate
# (see .cv_target_htilde() above); .retarget_boot_se() alone can't be reused
# for them since it hardcodes the plain H1/H0-only .final_target().
#
# Arguments:
#   Y, A     : scaled outcome and treatment
#   Q1W, Q0W : cross-fitted Q predictions
#   g1W      : propensity scores [length n]
#   H_tilde  : H_tilde covariate [length n]
#
# Returns: list(Q_star1, Q_star0)

.final_target_htilde <- function(Y, A, Q1W, Q0W, g1W, H_tilde) {
  g1W <- .clip_ps(g1W)

  d1f <- data.frame(Y = Y, H1 = A,     Ht = H_tilde, off1 = .logit(Q1W))
  d0f <- data.frame(Y = Y, H0 = 1 - A, Ht = H_tilde, off0 = .logit(Q0W))

  fit1f <- tryCatch(
    glm2(Y ~ -1 + H1 + Ht + offset(off1), family = quasibinomial,
         data = d1f, weights = 1 / g1W),
    error = function(e) NULL)
  fit0f <- tryCatch(
    glm2(Y ~ -1 + H0 + Ht + offset(off0), family = quasibinomial,
         data = d0f, weights = 1 / (1 - g1W)),
    error = function(e) NULL)

  if (is.null(fit1f) || is.null(fit0f)) return(NULL)

  eps1f  <- coef(fit1f)[1]; if (!is.finite(eps1f))  eps1f  <- 0
  eps1fg <- coef(fit1f)[2]; if (!is.finite(eps1fg)) eps1fg <- 0
  eps0f  <- coef(fit0f)[1]; if (!is.finite(eps0f))  eps0f  <- 0
  eps0fg <- coef(fit0f)[2]; if (!is.finite(eps0fg)) eps0fg <- 0

  list(
    Q_star1 = .expit(.logit(Q1W) + eps1f + eps1fg * H_tilde),
    Q_star0 = .expit(.logit(Q0W) + eps0f + eps0fg * H_tilde)
  )
}


# .retarget_boot_se_htilde()
#
# H_tilde-aware version of .retarget_boot_se(), for ctmle0()/ctmle0star().
# Same logic: Q, g, AND H_tilde all held fixed at their full-data values;
# only the targeting step (epsilon, epsilon_g) is refit per resample via
# .final_target_htilde(). See .retarget_boot_se() docstring for the general
# rationale and the fixed-nuisance caveat.
#
# Returns list(se, ci_lower, ci_upper, ate_boot, n_ok).

.retarget_boot_se_htilde <- function(Y, A, Q1W, Q0W, g1W, H_tilde, b_minus_a,
                                     alpha = 0.05, n_boot = 200) {
  n        <- length(Y)
  ate_boot <- rep(NA_real_, n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    tgt <- tryCatch(
      .final_target_htilde(Y[idx], A[idx], Q1W[idx], Q0W[idx],
                           g1W[idx], H_tilde[idx]),
      error = function(e) NULL)
    if (is.null(tgt)) next
    ate_b <- mean(tgt$Q_star1 - tgt$Q_star0)
    if (is.finite(ate_b)) ate_boot[b] <- ate_b
  }

  ate_boot <- ate_boot[is.finite(ate_boot)] * b_minus_a
  n_ok     <- length(ate_boot)
  if (n_ok < 2)
    return(list(se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                ate_boot = ate_boot, n_ok = n_ok))

  qs <- quantile(ate_boot, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  list(se = sd(ate_boot), ci_lower = qs[1], ci_upper = qs[2],
       ate_boot = ate_boot, n_ok = n_ok)
}


# .retarget_boot_se_wang()
#
# Retarget-only bootstrap for Wang-TMLE's parameterization, which is NOT the
# per-arm H1/H0 structure .retarget_boot_se() assumes: Wang-TMLE targets via
# a single clever covariate r = A - g(W) and offset regression
# Y ~ r - 1 + offset(Q0), giving one epsilon and beta1 = beta0 + epsilon
# (see wang_tmle() below) rather than separate Q_star1/Q_star0 arms.
#
# Q0, r, and beta0 are held FIXED at their full-data values (same fixed-
# nuisance convention as .retarget_boot_se()); only epsilon is refit per
# resample.
#
# Arguments:
#   Y, Q0, r  : scaled outcome, fixed outcome offset, fixed clever covariate
#   beta0     : fixed initial coefficient from the beta0 projection step
#   b_minus_a : Y-rescaling factor from .scale_y()
#
# Returns list(se, ci_lower, ci_upper, ate_boot, n_ok).

.retarget_boot_se_wang <- function(Y, Q0, r, beta0, b_minus_a,
                                   alpha = 0.05, n_boot = 200) {
  n        <- length(Y)
  ate_boot <- rep(NA_real_, n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    Yb  <- Y[idx]; rb <- r[idx]; Q0b <- Q0[idx]

    fit <- tryCatch(lm(Yb ~ rb - 1, offset = Q0b), error = function(e) NULL)
    if (is.null(fit) || anyNA(coef(fit))) next

    eps_b <- coef(fit)[1]
    if (is.finite(eps_b)) ate_boot[b] <- beta0 + eps_b
  }

  ate_boot <- ate_boot[is.finite(ate_boot)] * b_minus_a
  n_ok     <- length(ate_boot)
  if (n_ok < 2)
    return(list(se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                ate_boot = ate_boot, n_ok = n_ok))

  qs <- quantile(ate_boot, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  list(se = sd(ate_boot), ci_lower = qs[1], ci_upper = qs[2],
       ate_boot = ate_boot, n_ok = n_ok)
}


# compute_true_ate()

compute_true_ate <- function(model, snp_id, X, snp_ids) {
  snp_col <- match(snp_id, snp_ids)
  if (is.na(snp_col))
    stop("compute_true_ate: snp_id not found in snp_ids: ", snp_id)

  X1 <- X; X1[, snp_col] <- 1
  X0 <- X; X0[, snp_col] <- 0

  compute_mu <- function(Xmat) {
    mu <- numeric(nrow(Xmat))
    df <- as.data.frame(Xmat)
    colnames(df) <- snp_ids
    for (term in model) {
      vals <- df[, term$variants, drop = FALSE]
      mu   <- mu + term$coef * apply(vals, 1, prod)
    }
    mu
  }

  mean(compute_mu(X1) - compute_mu(X0))
}

# .retarget_boot_se()
#
# Fixed-nuisance ("re-target only") bootstrap SE/CI for a TMLE-type ATE.
#
# Q and g are held FIXED at their full-data fits; each resample re-runs ONLY
# the targeting step (epsilon) and re-averages. Captures variability of the
# targeting + averaging step, but by construction NOT the variability of
# estimating Q or g. Under strong LD the propensity model is the unstable
# piece (LASSO selection over near-collinear LD-block SNPs, near-positivity),
# so treat this as a CONTROL: if it stays near the EIF SE while a full-refit
# bootstrap matches the oracle SD, the missing variance is nuisance
# estimation, not targeting.
#
# Inputs are the already-computed length-n nuisance vectors from the point
# estimate (no refitting). b_minus_a rescales back to Y's original scale.
# NOTE: n_boot here is resamples PER replicate — distinct from the 200 MC
# replicates in ctmle_mc_analysis.R.
#
# Returns list(se, ci_lower, ci_upper, ate_boot, n_ok).

.retarget_boot_se <- function(Y, A, Q1W, Q0W, g1W, b_minus_a,
                              alpha = 0.05, n_boot = 500) {
  n        <- length(Y)
  ate_boot <- rep(NA_real_, n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)

    # Re-target on the resample using the FIXED nuisance predictions for
    # those rows. .final_target() fits epsilon per arm on the resampled
    # (Y, A, offset, weight) and returns Q_star for the resampled rows.
    # Duplicated rows are counted with multiplicity:  standard bootstrap.
    tgt <- tryCatch(
      .final_target(Y[idx], A[idx], Q1W[idx], Q0W[idx], g1W[idx]),
      error = function(e) NULL)
    if (is.null(tgt)) next

    ate_b <- mean(tgt$Q_star1 - tgt$Q_star0)
    if (is.finite(ate_b)) ate_boot[b] <- ate_b
  }

  ate_boot <- ate_boot[is.finite(ate_boot)] * b_minus_a
  n_ok     <- length(ate_boot)
  if (n_ok < 2)
    return(list(se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                ate_boot = ate_boot, n_ok = n_ok))

  qs <- quantile(ate_boot, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  list(se = sd(ate_boot), ci_lower = qs[1], ci_upper = qs[2],
       ate_boot = ate_boot, n_ok = n_ok)
}


# .g_refit_boot_se()
#
# Propensity-refit ("g-refit") bootstrap SE/CI for a TMLE-type ATE.
#
# Companion diagnostic to .retarget_boot_se(). That function holds Q AND g
# fixed and re-runs only the targeting step; if its SE still falls short of
# the oracle empirical SD, the natural next hypothesis (per the original
# LD-collinearity diagnosis: LASSO selecting arbitrarily among near-collinear
# SNPs) is that the missing variance lives in the propensity model SELECTION
# itself — i.e. which lambda / which SNPs cv.glmnet picks changes from one
# resample to the next under strong LD, and no fixed-g bootstrap can see that.
#
# This function re-runs the FULL cv.glmnet lambda search for g on each
# resample (capturing that selection instability), holds Q fixed at its
# original full-data cross-fit (Q/bias was never the diagnosed problem —
# no need to pay for a full outcome-model refit too), and re-targets with
# .final_target() using the freshly-fit g.
#
# Arguments:
#   Y, A, X   : scaled outcome, treatment, covariate matrix (pre-resampling)
#   Q1W, Q0W  : cross-fitted Q predictions from the ORIGINAL (non-bootstrap)
#               fit — held fixed across all resamples, same convention as
#               .retarget_boot_se()
#   b_minus_a : Y-rescaling factor from .scale_y()
#   n_boot    : bootstrap resamples PER replicate (same convention/caveat as
#               .retarget_boot_se() — distinct from the 200 MC replicates)
#
# Returns list(se, ci_lower, ci_upper, ate_boot, n_ok)

.g_refit_boot_se <- function(Y, A, X, Q1W, Q0W, b_minus_a,
                             alpha = 0.05, n_boot = 200) {
  n        <- length(Y)
  ate_boot <- rep(NA_real_, n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)

    # Full CV-lambda re-selection on the resample — this is the piece
    # .retarget_boot_se() cannot capture, since it never refits g.
    g_fit_b <- tryCatch(
      cv.glmnet(X[idx, , drop = FALSE], A[idx], family = "binomial", alpha = 1),
      error = function(e) NULL)
    if (is.null(g_fit_b)) next

    g1W_b <- tryCatch(
      .clip_ps(as.numeric(
        predict(g_fit_b, newx = X[idx, , drop = FALSE],
               s = "lambda.min", type = "response"))),
      error = function(e) NULL)
    if (is.null(g1W_b) || anyNA(g1W_b)) next

    # Q held FIXED at the original fit (indexed with multiplicity, as in
    # .retarget_boot_se()) — only g and epsilon vary across resamples.
    tgt_b <- tryCatch(
      .final_target(Y[idx], A[idx], Q1W[idx], Q0W[idx], g1W_b),
      error = function(e) NULL)
    if (is.null(tgt_b)) next

    ate_b <- mean(tgt_b$Q_star1 - tgt_b$Q_star0)
    if (is.finite(ate_b)) ate_boot[b] <- ate_b
  }

  ate_boot <- ate_boot[is.finite(ate_boot)] * b_minus_a
  n_ok     <- length(ate_boot)
  if (n_ok < 2)
    return(list(se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                ate_boot = ate_boot, n_ok = n_ok))

  qs <- quantile(ate_boot, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  list(se = sd(ate_boot), ci_lower = qs[1], ci_upper = qs[2],
       ate_boot = ate_boot, n_ok = n_ok)
}


# collect_results_snp() / collect_results()

collect_results_snp <- function(data, true_ate = NULL, alpha = 0.05) {

  # C-TMLE1 is pre-run once so its result can be reused without re-fitting
  ctmle1_res <- tryCatch(
    ctmle1(data, alpha = alpha),
    error = function(e) {
      message("Error in C-TMLE1: ", conditionMessage(e))
      .na_est("C-TMLE1", best_lambda = NULL)
    }
  )

  estimators <- list(
    "TMLE"            = function(d) tmle(d, alpha),
    "C-TMLE1"         = function(d) ctmle1_res,
    "C-TMLE0"         = function(d) ctmle0(d, alpha),
    "C-TMLE1-OAL"     = function(d) ctmle1_oal(d, alpha),
    "Scalable-CTMLE"  = function(d) scalable_ctmle(d, alpha),
    "Wang-TMLE"       = function(d) wang_tmle(d, alpha),
    "Wang-CTMLE"      = function(d) wang_ctmle(d, alpha)
  )

  rows <- lapply(names(estimators), function(nm) {
    res <- tryCatch(
      estimators[[nm]](data),
      error = function(e) {
        message("Error in ", nm, ": ", conditionMessage(e))
        .na_est(nm)
      }
    )
    bias <- if (!is.null(true_ate) && !is.na(res$est))
              res$est - true_ate else NA_real_
    data.frame(
      Algorithm      = nm,
      Estimate       = res$est,
      Standard.Error = res$se,
      CI_Lower       = res$ci_lower,
      CI_Upper       = res$ci_upper,
      true_ATE       = if (!is.null(true_ate)) true_ate else NA_real_,
      Bias           = bias,
      Bias.to.se     = if (!is.na(bias) && !is.na(res$se) && res$se > 0)
                         bias / res$se else NA_real_,
      pvalue         = res$pvalue,
      mean_IF        = if (!is.null(res$mean_IF)) res$mean_IF else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# Backward-compatible alias
collect_results <- function(data, true_ate = NULL, alpha = 0.05)
  collect_results_snp(data, true_ate = true_ate, alpha = alpha)


# MC_res_anchorgene()

MC_res_anchorgene <- function(combined_results) {
  has_true_ate <- "true_ATE" %in% colnames(combined_results) &&
                  !all(is.na(combined_results$true_ATE))

  summary_tbl <- combined_results |>
    dplyr::group_by(Algorithm) |>
    dplyr::summarize(
      N             = sum(!is.na(Estimate)),
      Mean_Est      = mean(Estimate,       na.rm = TRUE),
      Mean_Bias     = if (has_true_ate)
                        mean(Estimate - true_ATE, na.rm = TRUE)
                      else NA_real_,
      Mean_SE       = mean(Standard.Error, na.rm = TRUE),
      Sample_Var    = var(Estimate,        na.rm = TRUE),
      Mean_MSE      = if (has_true_ate)
                        mean((Estimate - true_ATE)^2, na.rm = TRUE)
                      else NA_real_,
      Coverage      = if (has_true_ate)
                        mean((true_ATE >= CI_Lower) &
                             (true_ATE <= CI_Upper), na.rm = TRUE)
                      else NA_real_,
      Mean_CI_Lower = mean(CI_Lower, na.rm = TRUE),
      Mean_CI_Upper = mean(CI_Upper, na.rm = TRUE),
      CI_2.5th      = quantile(Estimate, 0.025, na.rm = TRUE),
      CI_97.5th     = quantile(Estimate, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(., digits = 3)))

  if (!has_true_ate)
    summary_tbl <- dplyr::select(summary_tbl, -Mean_Bias, -Mean_MSE, -Coverage)

  summary_tbl
}


# Plotting helpers


.anchorgene_white_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      strip.text       = ggplot2::element_text(size = 11, face = "bold"),
      plot.title       = ggplot2::element_text(size = 13, face = "bold",
                                               hjust = 0.5)
    )
}

hist_per_est_anchorgene <- function(
    res,
    xlimits = range(res$Estimate, na.rm = TRUE),
    title   = "Distribution of ATE Estimates by Algorithm") {

  p <- ggplot2::ggplot(res, ggplot2::aes(x = Estimate, fill = Algorithm)) +
    ggplot2::geom_density(alpha = 0.4, color = "black", na.rm = TRUE) +
    ggplot2::scale_x_continuous(limits = xlimits) +
    ggplot2::facet_wrap(~Algorithm) +
    ggplot2::labs(title = title, x = "ATE Estimate", y = "Density") +
    .anchorgene_white_theme() +
    ggplot2::theme(legend.position = "none")

  if ("true_ATE" %in% colnames(res) && !all(is.na(res$true_ATE))) {
    mean_true <- res |>
      dplyr::group_by(Algorithm) |>
      dplyr::summarize(mean_true = mean(true_ATE, na.rm = TRUE), .groups = "drop")
    p <- p + ggplot2::geom_vline(
      data = mean_true,
      ggplot2::aes(xintercept = mean_true),
      color = "black", linetype = "dashed"
    )
  }
  p
}

hist_per_est_CI_anchorgene <- function(
    res,
    xlimits = c(min(res$CI_Lower, na.rm = TRUE),
                max(res$CI_Upper, na.rm = TRUE)),
    title   = "Distribution of CI Bounds by Algorithm") {

  res |>
    dplyr::select(Algorithm, CI_Lower, CI_Upper) |>
    tidyr::pivot_longer(cols = c(CI_Lower, CI_Upper),
                        names_to  = "CI_Bound",
                        values_to = "CI_Value") |>
    ggplot2::ggplot(ggplot2::aes(x = CI_Value,
                                 fill  = CI_Bound,
                                 color = CI_Bound)) +
    ggplot2::geom_density(alpha = 0.4, na.rm = TRUE) +
    ggplot2::scale_x_continuous(limits = xlimits) +
    ggplot2::facet_wrap(~Algorithm) +
    ggplot2::labs(title = title, x = "CI Bound Value", y = "Density",
                  fill = "CI Bound", color = "CI Bound") +
    .anchorgene_white_theme()
}

hist_all_anchorgene <- function(
    res,
    xlimits = range(res$Estimate, na.rm = TRUE),
    title   = "Distribution of ATE Estimates by Algorithm") {

  p <- ggplot2::ggplot(
    res,
    ggplot2::aes(x = Estimate, fill = Algorithm, color = Algorithm)
  ) +
    ggplot2::geom_density(alpha = 0.35, na.rm = TRUE) +
    ggplot2::scale_x_continuous(limits = xlimits) +
    ggplot2::labs(title = title, x = "ATE Estimate", y = "Density") +
    .anchorgene_white_theme()

  if ("true_ATE" %in% colnames(res) && !all(is.na(res$true_ATE))) {
    grand_mean_true <- mean(res$true_ATE, na.rm = TRUE)
    p <- p + ggplot2::geom_vline(
      xintercept = grand_mean_true,
      color = "black", linetype = "dashed", linewidth = 0.8
    )
  }
  p
}

plot_est_CI_anchorgene <- function(
    res,
    custom_title = "ATE Estimates with 95% CI",
    xlimits      = c(min(res$CI_Lower, na.rm = TRUE),
                     max(res$CI_Upper, na.rm = TRUE))) {

  if (any(duplicated(res$Algorithm))) {
    res <- res |>
      dplyr::group_by(Algorithm) |>
      dplyr::summarize(
        Estimate = mean(Estimate, na.rm = TRUE),
        CI_Lower = mean(CI_Lower, na.rm = TRUE),
        CI_Upper = mean(CI_Upper, na.rm = TRUE),
        .groups  = "drop"
      )
  }

  res <- res |>
    dplyr::mutate(CI_Width = CI_Upper - CI_Lower) |>
    dplyr::arrange(dplyr::desc(CI_Width))
  res$Algorithm <- factor(res$Algorithm, levels = unique(res$Algorithm))

  ggplot2::ggplot(res, ggplot2::aes(x = Algorithm, y = Estimate)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = CI_Lower, ymax = CI_Upper),
      width = 0.25, linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 3, color = "steelblue") +
    ggplot2::scale_y_continuous(limits = xlimits) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = custom_title, x = "", y = "ATE Estimate") +
    .anchorgene_white_theme() +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14)
    )
}