# estimators.R
#
# ESTIMATORS:
#   ctmle1()         — Collaborative LASSO, lambda by targeted CV
#   ctmle0()         — CV-selected lambda + extra clever covariate H~
#   ctmle0star()     — C-TMLE0 targeting with collaboratively-selected lambda
#   scalable_ctmle() — Pre-ordered partial-correlation forward selection
#   tmle()           — Standard TMLE
#   wang_tmle()      — Wang et al. TMLE
#   wang_ctmle()     — Wang et al. C-TMLE
#
# All shared helpers (package loading, .logit/.expit/.scale_y/.clip_ps,
# .W_matrix/.make_QX/.na_est, .crossfit_Q, .cv_target, .cv_target_htilde,
# .final_target, .collab_cv_loss, .ic_inference) live in utils.R.


source("src/utils.R")  # utils.R lives alongside this file under src/; every known
                        # caller (ctmle_mc_analysis.R, the training notebook, etc.)
                        # already runs with the working directory set to the repo
                        # root, so a plain relative path is simpler and more robust
                        # than introspecting the call stack. Do NOT add library() here.

# 1. C-TMLE1 — Collaborative LASSO (lambda selected by targeted CV)

ctmle1 <- function(data, alpha = 0.05, V = 10, cv_var = TRUE,
                   boot_var = TRUE, n_boot = 200) {

  sc <- .scale_y(data$Y)
  Y  <- sc$y_scaled
  A  <- data$A
  X  <- .W_matrix(data)
  n  <- nrow(data)

  if (length(unique(A)) < 2)
    return(.na_est("C-TMLE1", best_lambda = NA_real_))

  qx  <- .make_QX(A, X)
  cf  <- .crossfit_Q(qx$QX, qx$QX0, qx$QX1, Y, V = V)
  Q0W <- cf$Q0W
  Q1W <- cf$Q1W
  folds <- cf$folds

  ps_cv      <- cv.glmnet(X, A, family = "binomial", alpha = 1, nfolds = V)
  lambda_cv  <- ps_cv$lambda.min
  lambda_seq <- sort(
    ps_cv$glmnet.fit$lambda[ps_cv$glmnet.fit$lambda <= lambda_cv],
    decreasing = TRUE
  )
                                         
  best_lambda <- .collab_cv_loss(X, A, Y, Q1W, Q0W, folds, lambda_seq)

  # Cross-fitted propensity at the collaboratively selected lambda
  g1W <- .crossfit_g(X, A, folds, best_lambda)

  # Cross-fitted targeting epsilon (was: .final_target(), in-sample epsilon)
  tgt <- .cv_target(Y, A, Q1W, Q0W, g1W, folds) 
  if (is.null(tgt)) return(.na_est("C-TMLE1", best_lambda = best_lambda))

  ATE_n <- mean(tgt$Q_star1 - tgt$Q_star0)
  inf <- if (cv_var && boot_var) {
    boot <- .retarget_boot_se(Y, A, Q1W, Q0W, g1W, sc$b - sc$a,
                              alpha = alpha, n_boot = n_boot)
    ATE  <- ATE_n * (sc$b - sc$a)
    list(est = ATE, se = boot$se, ci_lower = boot$ci_lower,
         ci_upper = boot$ci_upper, pvalue = .pvalue(ATE, boot$se),
         mean_IF = NA_real_, n_boot_ok = boot$n_ok)
  } else if (cv_var) {
    .ic_inference_cv(A, Y, tgt$Q_star1, tgt$Q_star0, g1W, folds,
                     ATE_n, sc$b - sc$a, alpha)
  } else {
    .ic_inference(A, Y, tgt$Q_star1, tgt$Q_star0, g1W,
                  ATE_n, sc$b - sc$a, alpha)
  }
  c(inf, list(estimator = "C-TMLE1", best_lambda = best_lambda))
}

# 2. C-TMLE0 — CV-selected lambda + extra clever covariate H~

ctmle0 <- function(data, alpha = 0.05, V = 10, cv_var = TRUE,
                   boot_var = TRUE, n_boot = 200) {

  sc <- .scale_y(data$Y)
  Y  <- sc$y_scaled
  A  <- data$A
  X  <- .W_matrix(data)
  n  <- nrow(data)

  if (length(unique(A)) < 2) return(.na_est("C-TMLE0"))

  qx    <- .make_QX(A, X)
  cf    <- .crossfit_Q(qx$QX, qx$QX0, qx$QX1, Y, V = V)
  Q0W   <- cf$Q0W
  Q1W   <- cf$Q1W
  folds <- cf$folds

  ps_cv <- cv.glmnet(X, A, family = "binomial", alpha = 1, nfolds = V)
  g1W   <- .clip_ps(
    as.numeric(predict(ps_cv, newx = X, s = "lambda.min", type = "response"))
  )

  # H_tilde: gradient of propensity loss w.r.t. lambda
  lam_min <- ps_cv$lambda.min
  lam_d   <- min(exp(log(lam_min) + 1e-3), max(ps_cv$glmnet.fit$lambda))
  g1W_d   <- .clip_ps(
    as.numeric(predict(ps_cv, newx = X, s = lam_d, type = "response"))
  )
  dg      <- g1W_d - g1W
  H_tilde <- ((1 - A) / (1 - g1W)^2) * (-dg) + (A / g1W^2) * dg

  tgt <- if (ncol(X) == 1)
    .cv_target(Y, A, Q1W, Q0W, g1W, folds)
  else
    .cv_target_htilde(Y, A, Q1W, Q0W, g1W, H_tilde, folds)

  ATE_n <- mean(tgt$Q_star1 - tgt$Q_star0)
  inf   <- if (cv_var && boot_var) {
    # H_tilde is only defined/used when ncol(X) != 1 (see tgt above) — mirror
    # that fallback here so the bootstrap targets the same model the point
    # estimate actually used.
    boot <- if (ncol(X) == 1)
      .retarget_boot_se(Y, A, Q1W, Q0W, g1W, sc$b - sc$a,
                        alpha = alpha, n_boot = n_boot)
    else
      .retarget_boot_se_htilde(Y, A, Q1W, Q0W, g1W, H_tilde, sc$b - sc$a,
                               alpha = alpha, n_boot = n_boot)
    ATE <- ATE_n * (sc$b - sc$a)
    list(est = ATE, se = boot$se, ci_lower = boot$ci_lower,
         ci_upper = boot$ci_upper, pvalue = .pvalue(ATE, boot$se),
         mean_IF = NA_real_, n_boot_ok = boot$n_ok)
  } else if (cv_var) {
    .ic_inference_cv(A, Y, tgt$Q_star1, tgt$Q_star0, g1W, folds,
                     ATE_n, sc$b - sc$a, alpha)
  } else {
    .ic_inference(A, Y, tgt$Q_star1, tgt$Q_star0, g1W,
                  ATE_n, sc$b - sc$a, alpha)
  }
  c(inf, list(estimator = "C-TMLE0"))
}


# 3. C-TMLE0* — C-TMLE0 targeting with collaboratively-selected lambda

ctmle0star <- function(data, alpha = 0.05, V = 10, cv_var = TRUE,
                       boot_var = TRUE, n_boot = 200) {

  sc <- .scale_y(data$Y)
  Y  <- sc$y_scaled
  A  <- data$A
  X  <- .W_matrix(data)
  n  <- nrow(data)

  if (length(unique(A)) < 2) return(.na_est("C-TMLE0*"))

  qx    <- .make_QX(A, X)
  # Share folds across Q cross-fitting, collaborative CV, and final targeting
  folds <- caret::createFolds(Y, k = V, list = TRUE, returnTrain = FALSE)
  cf    <- .crossfit_Q(qx$QX, qx$QX0, qx$QX1, Y, V = V, folds = folds)
  Q0W   <- cf$Q0W
  Q1W   <- cf$Q1W

  ps_cv      <- cv.glmnet(X, A, family = "binomial", alpha = 1, nfolds = V)
  lambda_cv  <- ps_cv$lambda.min
  lambda_seq <- sort(
    ps_cv$glmnet.fit$lambda[ps_cv$glmnet.fit$lambda <= lambda_cv],
    decreasing = TRUE
  )

  best_lambda <- .collab_cv_loss(X, A, Y, Q1W, Q0W, folds, lambda_seq)
  if (!is.finite(best_lambda)) return(.na_est("C-TMLE0*"))

  g_model <- glmnet(X, A, family = "binomial", lambda = best_lambda)
  g1W     <- .clip_ps(as.numeric(predict(g_model, X, type = "response")))

  # H_tilde at best_lambda
  lam_d <- min(exp(log(best_lambda) + 1e-3), max(ps_cv$glmnet.fit$lambda))
  g1W_d <- .clip_ps(
    as.numeric(predict(ps_cv, newx = X, s = lam_d, type = "response"))
  )
  dg      <- g1W_d - g1W
  H_tilde <- ((1 - A) / (1 - g1W)^2) * (-dg) + (A / g1W^2) * dg

  tgt <- if (ncol(X) == 1)
    .cv_target(Y, A, Q1W, Q0W, g1W, folds)
  else
    .cv_target_htilde(Y, A, Q1W, Q0W, g1W, H_tilde, folds)

  ATE_n <- mean(tgt$Q_star1 - tgt$Q_star0)
  inf <- if (cv_var && boot_var) {
    boot <- if (ncol(X) == 1)
      .retarget_boot_se(Y, A, Q1W, Q0W, g1W, sc$b - sc$a,
                        alpha = alpha, n_boot = n_boot)
    else
      .retarget_boot_se_htilde(Y, A, Q1W, Q0W, g1W, H_tilde, sc$b - sc$a,
                               alpha = alpha, n_boot = n_boot)
    ATE <- ATE_n * (sc$b - sc$a)
    list(est = ATE, se = boot$se, ci_lower = boot$ci_lower,
         ci_upper = boot$ci_upper, pvalue = .pvalue(ATE, boot$se),
         mean_IF = NA_real_, n_boot_ok = boot$n_ok)
  } else if (cv_var) {
    .ic_inference_cv(A, Y, tgt$Q_star1, tgt$Q_star0, g1W, folds,
                     ATE_n, sc$b - sc$a, alpha)
  } else {
    .ic_inference(A, Y, tgt$Q_star1, tgt$Q_star0, g1W,
                  ATE_n, sc$b - sc$a, alpha)
  }
  c(inf, list(estimator = "C-TMLE0*"))
}


# 4. Scalable C-TMLE — partial-correlation pre-ordering, O(p) selection

scalable_ctmle <- function(data, alpha = 0.05, V = 10, cv_var = TRUE,
                           boot_var = TRUE, n_boot = 200) {

  sc <- .scale_y(data$Y)
  Y  <- sc$y_scaled
  A  <- data$A
  X  <- .W_matrix(data)
  n  <- nrow(data)

  if (ncol(X) < 2) {
    res <- tmle(data, alpha = alpha)
    res$estimator <- "Scalable-CTMLE"
    return(res)
  }
  if (length(unique(A)) < 2)
    return(.na_est("Scalable-CTMLE"))

  qx    <- .make_QX(A, X)
  cf    <- .crossfit_Q(qx$QX, qx$QX0, qx$QX1, Y, V = V)
  Q0W   <- cf$Q0W
  Q1W   <- cf$Q1W
  Q_cur <- cf$Q_full
  folds <- cf$folds                        # reuse split for g + epsilon

  residual  <- Y - Q_cur
  cor_order <- order(abs(cor(X, residual)), decreasing = TRUE)
  X_ordered <- X[, cor_order, drop = FALSE]

  best_k    <- NA_integer_
  best_loss <- Inf

  # Model SELECTION unchanged: in-sample g here is fine, it only ranks models.
  for (k in seq(2, min(ncol(X_ordered), 50))) {
    Xk    <- X_ordered[, seq_len(k), drop = FALSE]
    g_fit <- tryCatch(
      cv.glmnet(Xk, A, family = "binomial", alpha = 1, nfolds = V),
      error = function(e) NULL)
    if (is.null(g_fit)) next

    g1W_k <- .clip_ps(
      as.numeric(predict(g_fit, newx = Xk, s = "lambda.min", type = "response"))
    )

    d1 <- data.frame(Y = Y, H1 = A,     off = .logit(Q1W))
    d0 <- data.frame(Y = Y, H0 = 1 - A, off = .logit(Q0W))
    fit1 <- tryCatch(
      glm2(Y ~ -1 + H1 + offset(off), family = quasibinomial,
           data = d1, weights = 1 / g1W_k),
      error = function(e) NULL)
    fit0 <- tryCatch(
      glm2(Y ~ -1 + H0 + offset(off), family = quasibinomial,
           data = d0, weights = 1 / (1 - g1W_k)),
      error = function(e) NULL)
    if (is.null(fit1) || is.null(fit0)) next

    eps1 <- coef(fit1)[1]; eps0 <- coef(fit0)[1]
    Qs1  <- .expit(.logit(Q1W) + eps1)
    Qs0  <- .expit(.logit(Q0W) + eps0)
    loss <- mean((Y - (A * Qs1 + (1 - A) * Qs0))^2)

    if (loss < best_loss) { best_loss <- loss; best_k <- k }
  }

  if (is.na(best_k)) return(.na_est("Scalable-CTMLE"))   # degenerate: no model fit

  # Winning model = top-best_k residual-ordered covariates.
  Xbest   <- X_ordered[, seq_len(best_k), drop = FALSE]
  ps_best <- cv.glmnet(Xbest, A, family = "binomial", alpha = 1, nfolds = V)

  if (cv_var) {
    g1W <- .crossfit_g(Xbest, A, folds, ps_best$lambda.min)     # OUT-OF-FOLD g
    tgt <- .cv_target(Y, A, Q1W, Q0W, g1W, folds)               # OUT-OF-FOLD eps
    if (is.null(tgt)) return(.na_est("Scalable-CTMLE"))
    ATE_n <- mean(tgt$Q_star1 - tgt$Q_star0)

    if (boot_var) {
      boot <- .retarget_boot_se(Y, A, Q1W, Q0W, g1W, sc$b - sc$a,
                                alpha = alpha, n_boot = n_boot)
      ATE  <- ATE_n * (sc$b - sc$a)
      inf  <- list(est = ATE, se = boot$se, ci_lower = boot$ci_lower,
                  ci_upper = boot$ci_upper, pvalue = .pvalue(ATE, boot$se),
                  mean_IF = NA_real_, n_boot_ok = boot$n_ok)
    } else {
      inf <- .ic_inference_cv(A, Y, tgt$Q_star1, tgt$Q_star0, g1W, folds,
                              ATE_n, sc$b - sc$a, alpha)
    }
  } else {
    g1W <- .clip_ps(as.numeric(
             predict(ps_best, newx = Xbest, s = "lambda.min", type = "response")))
    tgt <- .final_target(Y, A, Q1W, Q0W, g1W)
    if (is.null(tgt)) return(.na_est("Scalable-CTMLE"))
    ATE_n <- mean(tgt$Q_star1 - tgt$Q_star0)
    inf   <- .ic_inference(A, Y, tgt$Q_star1, tgt$Q_star0, g1W,
                           ATE_n, sc$b - sc$a, alpha)
  }

  c(inf, list(estimator = "Scalable-CTMLE"))
}

# 5. Standard TMLE  (+ CV-Var inference to match the collaborative estimators)
#
# cv_var = TRUE  -> cross-fit g (.crossfit_g) + cross-fit epsilon (.cv_target)
#                   + CV-Var SE (.ic_inference_cv). This is the coverage fix.
# cv_var = FALSE -> original in-sample g + .final_target + plug-in EIF SE.
# The boot/.retarget_boot_se path is kept only as a diagnostic control.

tmle <- function(data, alpha = 0.05, V = 10, cv_var = TRUE,
                 boot_var = TRUE, n_boot = 200) {

  sc <- .scale_y(data$Y)
  Y  <- sc$y_scaled
  A  <- data$A
  X  <- .W_matrix(data)

  qx  <- .make_QX(A, X)
  cf  <- .crossfit_Q(qx$QX, qx$QX0, qx$QX1, Y, V = V)
  Q0W <- cf$Q0W
  Q1W <- cf$Q1W
  folds <- cf$folds                       # reuse the SAME split for g and epsilon

  if (cv_var) {
    # lambda chosen by CV, but g predicted OUT-OF-FOLD at that lambda
    ps_cv <- cv.glmnet(X, A, family = "binomial", alpha = 1, nfolds = V)
    g1W   <- .crossfit_g(X, A, folds, ps_cv$lambda.min)
    tgt   <- .cv_target(Y, A, Q1W, Q0W, g1W, folds)
    if (is.null(tgt)) return(.na_est("TMLE"))
    ATE_n <- mean(tgt$Q_star1 - tgt$Q_star0)

    if (boot_var) {
      # PRODUCTION FIX (pilot-validated): CV-Var's plug-in SE undercovers
      # under strong LD (48% coverage on rs0254) because it treats g as
      # fixed. Retarget-only bootstrap holds Q1W/Q0W/g1W fixed at THIS
      # fit's values and re-runs only the targeting step across resamples,
      # which recovered 84% coverage in pilot_targeted_bootstrap.R — still
      # short of the oracle SD, but the best validated option cheap enough
      # for production (no CV refitting inside the bootstrap).
      boot <- .retarget_boot_se(Y, A, Q1W, Q0W, g1W, sc$b - sc$a,
                                alpha = alpha, n_boot = n_boot)
      ATE  <- ATE_n * (sc$b - sc$a)
      inf  <- list(
        est       = ATE,
        se        = boot$se,
        ci_lower  = boot$ci_lower,
        ci_upper  = boot$ci_upper,
        pvalue    = .pvalue(ATE, boot$se),
        mean_IF   = NA_real_,
        n_boot_ok = boot$n_ok
      )
    } else {
      inf <- .ic_inference_cv(A, Y, tgt$Q_star1, tgt$Q_star0, g1W, folds,
                              ATE_n, sc$b - sc$a, alpha)
    }
  } else {
    ps_cv <- cv.glmnet(X, A, family = "binomial", alpha = 1, nfolds = V)
    g1W   <- .clip_ps(as.numeric(
               predict(ps_cv, newx = X, s = "lambda.min", type = "response")))
    tgt   <- .final_target(Y, A, Q1W, Q0W, g1W)
    if (is.null(tgt)) return(.na_est("TMLE"))
    ATE_n <- mean(tgt$Q_star1 - tgt$Q_star0)
    inf   <- .ic_inference(A, Y, tgt$Q_star1, tgt$Q_star0, g1W,
                           ATE_n, sc$b - sc$a, alpha)
  }

  c(inf, list(estimator = "TMLE"))
}

# 6. Wang-TMLE

wang_tmle <- function(data, alpha = 0.05, V = 10, cv_var = TRUE,
                      boot_var = TRUE, n_boot = 200) {

  sc  <- .scale_y(data$Y)
  Y   <- sc$y_scaled
  A   <- data$A
  X   <- .W_matrix(data)
  n   <- length(Y)
  bma <- sc$b - sc$a

  if (length(unique(A)) < 2) return(.na_est("Wang-TMLE"))

  # 1. Cross-fitted background Q (prevents overfitting absorption of A)
  QX_full  <- model.matrix(~ A + ., data = data.frame(A = A, X))
  QX0_full <- model.matrix(~ A + ., data = data.frame(A = 0, X = X))
  folds    <- caret::createFolds(Y, k = V, list = TRUE, returnTrain = FALSE)
  QA0      <- numeric(n)
  Q0       <- numeric(n)

  for (v in seq_len(V)) {
    tr  <- setdiff(seq_len(n), folds[[v]])
    val <- folds[[v]]
    Qf  <- tryCatch(
      cv.glmnet(QX_full[tr, ], Y[tr], alpha = 1, nfolds = V),
      error = function(e) NULL)
    if (is.null(Qf)) next
    pred_QA0 <- tryCatch(
      as.numeric(predict(Qf, newx = QX0_full[val, ], s = "lambda.min")),
      error = function(e) NULL)
    pred_Q0  <- tryCatch(
      as.numeric(predict(Qf, newx = QX_full[val, ],  s = "lambda.min")),
      error = function(e) NULL)
    if (!is.null(pred_QA0) && !anyNA(pred_QA0)) QA0[val] <- pred_QA0
    if (!is.null(pred_Q0)  && !anyNA(pred_Q0))  Q0[val]  <- pred_Q0
  }
  
  # Check if cross-fitting produced any valid predictions
  if (all(QA0 == 0) || all(Q0 == 0)) return(.na_est("Wang-TMLE"))

  # Targeting offset Q0.
  #   cv_var = TRUE  -> keep the CROSS-FIT Q0 from the loop above. The old
  #                     in-sample refit re-introduced the outcome-model overfit
  #                     that shrinks (Y - Q_star) and deflates IC under strong
  #                     LD -- the same mechanism that broke the plug-in SE.
  #   cv_var = FALSE -> original behaviour: overwrite Q0 with a full-data
  #                     in-sample fit.
  if (!cv_var) {
    Q_init <- tryCatch(
      cv.glmnet(QX_full, Y, alpha = 1),
      error = function(e) NULL)
    if (is.null(Q_init)) return(.na_est("Wang-TMLE"))

    Q0_new <- tryCatch(
      as.numeric(predict(Q_init, newx = QX_full, s = "lambda.min")),
      error = function(e) NULL)
    if (is.null(Q0_new) || anyNA(Q0_new)) return(.na_est("Wang-TMLE"))
    Q0 <- Q0_new
  }

  beta0_fit <- tryCatch(
    lm(Y ~ A, offset = QA0),
    error = function(e) NULL)
  if (is.null(beta0_fit) || anyNA(coef(beta0_fit))) return(.na_est("Wang-TMLE"))
  beta0  <- coef(beta0_fit)[["A"]]

  # 2. Propensity g_n(W) = E[A | W]
  #   cv_var = TRUE  -> OUT-OF-FOLD g over the SAME folds used for Q, so the
  #                     propensity for an observation never comes from a model
  #                     that has already seen it. This is the piece that keeps
  #                     r = A - g (and hence IC) from collapsing under strong LD.
  #   cv_var = FALSE -> original in-sample cv.glmnet fit/predict.
  gAW <- if (ncol(X) == 0) {
    rep(mean(A), n)
  } else if (ncol(X) == 1) {
    lm_g <- tryCatch(lm(A ~ X), error = function(e) NULL)
    if (is.null(lm_g)) return(.na_est("Wang-TMLE"))
    fitted(lm_g)
  } else if (cv_var) {
    gcf <- numeric(n)
    for (v in seq_len(V)) {
      tr  <- setdiff(seq_len(n), folds[[v]])
      val <- folds[[v]]
      gf  <- tryCatch(cv.glmnet(X[tr, , drop = FALSE], A[tr], alpha = 1),
                      error = function(e) NULL)
      pg  <- if (is.null(gf)) NULL else tryCatch(
        as.numeric(predict(gf, newx = X[val, , drop = FALSE], s = "lambda.min")),
        error = function(e) NULL)
      gcf[val] <- if (is.null(pg) || anyNA(pg)) mean(A[tr]) else pg
    }
    gcf
  } else {
    g_fit <- tryCatch(cv.glmnet(X, A, alpha = 1), error = function(e) NULL)
    if (is.null(g_fit)) return(.na_est("Wang-TMLE"))
    pred_g <- tryCatch(
      as.numeric(predict(g_fit, newx = X, s = "lambda.min")),
      error = function(e) NULL)
    if (is.null(pred_g) || anyNA(pred_g)) return(.na_est("Wang-TMLE"))
    pred_g
  }
  if (anyNA(gAW)) return(.na_est("Wang-TMLE"))
  gAW <- .clip_ps(gAW)

  # 3. Clever covariate r = A - g_n(W)
  r <- A - gAW

  # 4. epsilon-regression
  eps_fit  <- tryCatch(lm(Y ~ r - 1, offset = Q0), error = function(e) NULL)
  if (is.null(eps_fit) || anyNA(coef(eps_fit))) return(.na_est("Wang-TMLE"))
  
  eps_coef <- summary(eps_fit)$coef
  if (is.null(eps_coef) || nrow(eps_coef) == 0) return(.na_est("Wang-TMLE"))
  eps_n    <- eps_coef[1, 1]  # Estimate column
  eps_p    <- eps_coef[1, 4]  # p-value column

  # 5. Update
  Q_star <- tryCatch(predict(eps_fit), error = function(e) NULL)
  if (is.null(Q_star) || anyNA(Q_star)) return(.na_est("Wang-TMLE"))
  beta1  <- beta0 + eps_n

  # 6. Influence curve and inference
  dDh     <- sum(A * r) / n
  dDh_min <- var(A) / sqrt(n)
  if (!is.finite(dDh) || abs(dDh) < dDh_min) return(.na_est("Wang-TMLE"))

  IC <- r * (Y - Q_star) / dDh
  if (any(is.na(IC))) return(.na_est("Wang-TMLE"))

  IC_var <- mean(IC^2)
  RSS    <- sum((Q_star - Y)^2)

  est <- beta1 * bma

  if (boot_var) {
    # PRODUCTION FIX, Wang-TMLE variant: r, Q0, beta0 held fixed at this
    # fit's values; only epsilon is refit per resample. Wang-TMLE's
    # parameterization (single offset regression, no per-arm Q_star1/
    # Q_star0) means it can't reuse .retarget_boot_se() directly — see
    # .retarget_boot_se_wang() in utils.R.
    boot <- .retarget_boot_se_wang(Y, Q0, r, beta0, bma,
                                   alpha = alpha, n_boot = n_boot)
    SE   <- boot$se
    list(
      est       = est,
      se        = SE,
      ci_lower  = boot$ci_lower,
      ci_upper  = boot$ci_upper,
      pvalue    = .pvalue(est, SE),
      estimator = "Wang-TMLE",
      IC_var    = IC_var,
      RSS       = RSS,
      eps       = eps_n,
      eps_p     = eps_p,
      n_boot_ok = boot$n_ok
    )
  } else {
    SE_sc <- sqrt(IC_var / n)
    SE    <- SE_sc * bma
    z     <- qnorm(1 - alpha / 2)

    list(
      est       = est,
      se        = SE,
      ci_lower  = est - z * SE,
      ci_upper  = est + z * SE,
      pvalue    = 2 * pnorm(-abs(est / SE)),
      estimator = "Wang-TMLE",
      IC_var    = IC_var,
      RSS       = RSS,
      eps       = eps_n,
      eps_p     = eps_p
    )
  }
}


# 7. Wang-CTMLE

wang_ctmle <- function(data, alpha = 0.05, cor.cut = 0.7, V = 10,
                       boot_var = TRUE, n_boot = 200) {

  sc  <- .scale_y(data$Y)
  Y   <- sc$y_scaled
  A   <- data$A
  X   <- .W_matrix(data)
  n   <- length(Y)
  p   <- ncol(X)
  bma <- sc$b - sc$a

  if (length(unique(A)) < 2 || p == 0) return(.na_est("Wang-CTMLE"))

  K.max    <- min(p, 10L)
  QX_full  <- model.matrix(~ A + ., data = data.frame(A = A, X))
  QX0_full <- model.matrix(~ A + ., data = data.frame(A = 0, X = X))

  # Cross-fitted background Q + beta0 (prevents in-sample overfit absorbing A)
  folds_cf <- caret::createFolds(Y, k = V, list = TRUE, returnTrain = FALSE)
  QA0_cf   <- numeric(n)
  Q0n_cf   <- numeric(n)

  for (v in seq_len(V)) {
    tr  <- setdiff(seq_len(n), folds_cf[[v]])
    val <- folds_cf[[v]]
    Qf  <- tryCatch(
      cv.glmnet(QX_full[tr, ], Y[tr], alpha = 1, nfolds = V),
      error = function(e) NULL)
    if (is.null(Qf)) next
    pred_QA0_cf <- tryCatch(
      as.numeric(predict(Qf, newx = QX0_full[val, ], s = "lambda.min")),
      error = function(e) NULL)
    pred_Q0_cf  <- tryCatch(
      as.numeric(predict(Qf, newx = QX_full[val, ],  s = "lambda.min")),
      error = function(e) NULL)
    if (!is.null(pred_QA0_cf) && !anyNA(pred_QA0_cf)) QA0_cf[val] <- pred_QA0_cf
    if (!is.null(pred_Q0_cf)  && !anyNA(pred_Q0_cf))  Q0n_cf[val] <- pred_Q0_cf
  }

  if (all(QA0_cf == 0) || all(Q0n_cf == 0)) return(.na_est("Wang-CTMLE"))
  Q0n <- Q0n_cf

  beta0_fit <- tryCatch(lm(Y ~ A, offset = QA0_cf), error = function(e) NULL)
  if (is.null(beta0_fit) || anyNA(coef(beta0_fit))) return(.na_est("Wang-CTMLE"))
  beta0 <- coef(beta0_fit)[["A"]]

  # Inner TMLE step (Wang parameterisation: single clever covariate r = A - g)
  .run_tmle <- function(gW_in, Q0_in, beta0_in, Y_in, A_in) {
    n_in <- length(Y_in)
    gW   <- .clip_ps(gW_in)
    r    <- A_in - gW
    fit  <- tryCatch(lm(Y_in ~ r - 1, offset = Q0_in), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    eps <- coef(fit)[["r"]]
    if (!is.finite(eps)) return(NULL)
    Q1n   <- Q0_in + eps * r
    beta1 <- beta0_in + eps
    mse   <- mean((Y_in - Q1n)^2)
    dDh   <- mean(A_in * r)
    if (!is.finite(dDh) || abs(dDh) < 1e-10) return(NULL)
    IC     <- r * (Y_in - Q1n) / dDh
    sigma2 <- mean(IC^2) / n_in
    list(Q1n = Q1n, gW = gW, beta = beta1,
         mse = mse, sigma2 = sigma2, pMSE = mse + sigma2, IC = IC)
  }

  # Forward selection of clever covariates (correlation-pruned, pMSE-greedy)
  gW_intercept <- .clip_ps(rep(mean(A), n))
  tmle0 <- .run_tmle(gW_intercept, Q0n, beta0, Y, A)
  if (is.null(tmle0)) return(.na_est("Wang-CTMLE"))

  candidates      <- vector("list", K.max + 1L)
  candidates[[1]] <- list(
    gW = gW_intercept, Q0n = Q0n, beta0 = beta0,
    pMSE = tmle0$pMSE, tmle_out = tmle0,
    vars = character(0), n_clever = 1L
  )

  current_Q0n <- Q0n
  current_b0  <- beta0
  available   <- seq_len(p)
  in_model    <- integer(0)
  n_clever    <- 1L

  for (k in seq_len(K.max)) {
    prev_cand <- candidates[[k]]
    prev_pMSE <- prev_cand$pMSE

    # Correlation-based pruning of available covariates
    if (length(in_model) > 0 && length(available) > 0) {
      cors    <- abs(cor(X[, available, drop = FALSE],
                         X[, in_model,  drop = FALSE]))
      max_cor <- if (length(in_model) == 1)
                   as.numeric(cors)
                 else
                   apply(cors, 1, max)
      available <- available[max_cor <= cor.cut]
    }

    if (length(available) == 0) {
      candidates[[k + 1L]] <- prev_cand
      next
    }

    best_pMSE <- Inf; best_idx <- NA_integer_; best_tmle <- NULL

    for (j in available) {
      Xg   <- X[, c(in_model, j), drop = FALSE]
      gf   <- tryCatch(lm(A ~ Xg), error = function(e) NULL)
      if (is.null(gf)) next
      t_j  <- .run_tmle(fitted(gf), current_Q0n, current_b0, Y, A)
      if (is.null(t_j)) next
      if (t_j$pMSE < best_pMSE) {
        best_pMSE <- t_j$pMSE; best_idx <- j; best_tmle <- t_j
      }
    }

    if (!is.na(best_idx) && best_pMSE < prev_pMSE) {
      in_model  <- c(in_model, best_idx)
      available <- setdiff(available, best_idx)
      gW_new    <- fitted(lm(A ~ X[, in_model, drop = FALSE]))
      candidates[[k + 1L]] <- list(
        gW = gW_new, Q0n = current_Q0n, beta0 = current_b0,
        pMSE = best_pMSE, tmle_out = best_tmle,
        vars = colnames(X)[in_model], n_clever = n_clever
      )
    } else {
      n_clever    <- n_clever + 1L
      current_Q0n <- prev_cand$tmle_out$Q1n
      current_b0  <- prev_cand$tmle_out$beta
      t_retry     <- .run_tmle(prev_cand$gW, current_Q0n, current_b0, Y, A)
      candidates[[k + 1L]] <- if (!is.null(t_retry))
        list(gW = prev_cand$gW, Q0n = current_Q0n, beta0 = current_b0,
             pMSE = t_retry$pMSE, tmle_out = t_retry,
             vars = prev_cand$vars, n_clever = n_clever)
      else
        prev_cand
    }
  }

  # CV model selection over candidates. The Q model is cached ONCE per fold
  # and reused across every candidate (only g varies across candidates).
  folds     <- caret::createFolds(A, k = V, list = TRUE, returnTrain = FALSE)
  best_star <- Inf
  best_k    <- 1L

  fold_Q_cache <- vector("list", V)

  for (v in seq_len(V)) {
    tr  <- setdiff(seq_len(n), folds[[v]])
    val <- folds[[v]]

    QX_tr   <- model.matrix(~ A + ., data = data.frame(A = A[tr],  X[tr,  , drop = FALSE]))
    QX0_val <- model.matrix(~ A + ., data = data.frame(A = 0,      X[val, , drop = FALSE]))
    QX_val  <- model.matrix(~ A + ., data = data.frame(A = A[val], X[val, , drop = FALSE]))

    Qf_v <- tryCatch(
      cv.glmnet(QX_tr, Y[tr], alpha = 1, nfolds = V),
      error = function(e) NULL)
    if (is.null(Qf_v)) next

    pred_QA0_val <- tryCatch(
      as.numeric(predict(Qf_v, newx = QX0_val, s = "lambda.min")),
      error = function(e) NULL)
    pred_Q0_val  <- tryCatch(
      as.numeric(predict(Qf_v, newx = QX_val,  s = "lambda.min")),
      error = function(e) NULL)
    if (is.null(pred_QA0_val) || anyNA(pred_QA0_val)) next
    if (is.null(pred_Q0_val)  || anyNA(pred_Q0_val))  next

    b0_fit_val <- tryCatch(
      lm(Yv ~ Av, data = data.frame(Yv = Y[val], Av = A[val]), offset = pred_QA0_val),
      error = function(e) NULL)
    if (is.null(b0_fit_val)) next

    b0_val <- tryCatch(coef(b0_fit_val)[["Av"]], error = function(e) NA_real_)
    if (!is.finite(b0_val)) next

    fold_Q_cache[[v]] <- list(tr = tr, val = val, Q0_val = pred_Q0_val, b0_val = b0_val)
  }

  for (k in seq_along(candidates)) {
    cand <- candidates[[k]]
    if (is.null(cand)) next

    cv_mse <- 0; n_valid <- 0

    for (v in seq_len(V)) {
      fq <- fold_Q_cache[[v]]
      if (is.null(fq)) next
      tr     <- fq$tr
      val    <- fq$val
      Q0_val <- fq$Q0_val
      b0_val <- fq$b0_val

      g_vars <- cand$vars
      gW_val <- if (length(g_vars) == 0) {
        rep(mean(A[tr]), length(val))
      } else {
        gf_v <- tryCatch(
          lm(A[tr] ~ X[tr, g_vars, drop = FALSE]),
          error = function(e) NULL)
        if (is.null(gf_v)) next
        as.numeric(cbind(1, X[val, g_vars, drop = FALSE]) %*% coef(gf_v))
      }

      t_v <- .run_tmle(gW_val, Q0_val, b0_val, Y[val], A[val])
      if (is.null(t_v)) next
      cv_mse  <- cv_mse + t_v$mse * length(val)
      n_valid <- n_valid + length(val)
    }

    if (n_valid == 0) next
    cv_mse    <- cv_mse / n_valid
    size_k    <- length(cand$vars)
    pMSE_star <- n * log(cv_mse + cand$tmle_out$sigma2) + size_k * log(n)

    if (pMSE_star < best_star) {
      best_star <- pMSE_star
      best_k    <- k
    }
  }

  # Final cross-fit g at the winning variable set, then retarget
  best_vars     <- candidates[[best_k]]$vars
  gW_cf         <- .crossfit_g_lm(X, A, folds_cf, best_vars)
  beta0_best_in <- candidates[[best_k]]$beta0
  Q0_best_in    <- candidates[[best_k]]$Q0n

  tmle_best <- .run_tmle(gW_cf, Q0_best_in, beta0_best_in, Y, A)
  if (is.null(tmle_best))
    # Cross-fit retargeting failed (e.g. degenerate dDh) -- fall back to the
    # in-sample candidate rather than erroring the whole estimator out.
    tmle_best <- candidates[[best_k]]$tmle_out
  IC  <- tmle_best$IC
  ATE <- tmle_best$beta * bma

  if (boot_var) {
    r_best <- A - gW_cf
    boot <- .retarget_boot_se_wang(Y, Q0_best_in, r_best, beta0_best_in, bma,
                                   alpha = alpha, n_boot = n_boot)
    list(
      est          = ATE,
      se           = boot$se,
      ci_lower     = boot$ci_lower,
      ci_upper     = boot$ci_upper,
      pvalue       = .pvalue(ATE, boot$se),
      estimator    = "Wang-CTMLE",
      n_boot_ok    = boot$n_ok,
      best_k       = best_k,
      n_vars_final = length(best_vars)
    )
  } else {
    SE <- sqrt(mean(IC^2) / n) * bma
    z  <- qnorm(1 - alpha / 2)

    list(
      est          = ATE,
      se           = SE,
      ci_lower     = ATE - z * SE,
      ci_upper     = ATE + z * SE,
      pvalue       = .pvalue(ATE, SE),
      estimator    = "Wang-CTMLE",
      best_k       = best_k,
      n_vars_final = length(best_vars)
    )
  }
}


# 8. C-TMLE1-OAL — Collaborative-controlled, cross-fit outcome-adaptive LASSO
#
# Intended to replicate Wyss et al. (2024, AJE), Model 8 ("Collaborative-
# Controlled CF OAL"): an outcome LASSO selects confounders (penalty.factor
# = 0 for treatment, line below), those variables are forced unpenalized
# into an adaptive LASSO for treatment, the degree of undersmoothing is
# chosen by collaborative-controlled targeted learning rather than CV
# prediction error, and the final propensity score is cross-fitted.


ctmle1_oal <- function(data, alpha = 0.05, V = 10, cv_var = TRUE,
                       boot_var = TRUE, n_boot = 200) {

  sc <- .scale_y(data$Y)
  Y  <- sc$y_scaled
  A  <- data$A
  X  <- .W_matrix(data)
  n  <- nrow(data)
  p  <- ncol(X)

  if (length(unique(A)) < 2)
    return(.na_est("C-TMLE1-OAL",
                   best_lambda     = NA_real_,
                   n_forced        = NA_integer_,
                   chosen_vars_out = character(0)))

  qx  <- .make_QX(A, X)
  cf  <- .crossfit_Q(qx$QX, qx$QX0, qx$QX1, Y, V = V)
  Q0W <- cf$Q0W
  Q1W <- cf$Q1W
  folds <- cf$folds                       # CHANGED: one split for Q, g, targeting, variance

  # Outcome LASSO — variable selection for the OAL penalty vector
  # A is forced in (penalty.factor = 0 for its column) so selected W columns
  # reflect outcome associations not absorbed by the treatment effect.
  # model.matrix layout: (Intercept), A, W1, W2, ...
  # penalty.factor covers non-intercept columns only.
  pf_outcome    <- rep(1, ncol(qx$QX) - 1)   # drop intercept column
  pf_outcome[1] <- 0                           # first non-intercept col = A

  out_lasso <- tryCatch(
    cv.glmnet(
      x              = qx$QX[, -1, drop = FALSE],
      y              = Y,
      alpha          = 1,
      nfolds         = 5,
      penalty.factor = pf_outcome,
      standardize    = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(out_lasso))
    return(.na_est("C-TMLE1-OAL",
                   best_lambda     = NA_real_,
                   n_forced        = NA_integer_,
                   chosen_vars_out = character(0)))

  out_coef      <- as.numeric(coef(out_lasso, s = "lambda.min"))
  coef_names    <- rownames(coef(out_lasso, s = "lambda.min"))
  w_col_names   <- colnames(X)
  coef_w_idx    <- which(coef_names %in% w_col_names)
  chosen_vars_out <- coef_names[coef_w_idx][out_coef[coef_w_idx] != 0]

  pf_treatment <- rep(1, p)
  names(pf_treatment) <- colnames(X)
  pf_treatment[colnames(X) %in% chosen_vars_out] <- 0

  oal_fit <- tryCatch(
    cv.glmnet(
      x              = X,
      y              = A,
      family         = "binomial",
      alpha          = 1,
      nfolds         = V,
      penalty.factor = pf_treatment,
      standardize    = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(oal_fit))
    return(.na_est("C-TMLE1-OAL",
                   best_lambda     = NA_real_,
                   n_forced        = length(chosen_vars_out),
                   chosen_vars_out = chosen_vars_out))
                                          # CHANGED: removed separate caret::createFolds(A, ...)
  oal_ps_fn <- function(X_tr, A_tr, lam)
    glmnet(
      x              = X_tr,
      y              = A_tr,
      family         = "binomial",
      lambda         = lam,
      penalty.factor = pf_treatment,
      standardize    = TRUE
    )

  # Undersmoothing search direction (bug fix #1)
  # Wyss et al. (2024, Appendix S1) start at lambda.min and walk toward
  # SMALLER lambda (less regularization), stopping once the same-sample AUC
  # of the treatment model exceeds 0.999. path_lambda is in glmnet's default
  # decreasing order, so "less regularization" means higher column index.
  path_lambda <- oal_fit$glmnet.fit$lambda
  lambda_cv   <- oal_fit$lambda.min
  pos_start   <- which.min(abs(path_lambda - lambda_cv))

  # Same-sample predicted probabilities across the whole path — cheap, reuses
  # the glmnet object already fit inside cv.glmnet (no refitting needed).
  gns_path <- predict(oal_fit$glmnet.fit, newx = X, type = "response")
  auc_path <- vapply(seq_along(path_lambda),
                     function(j) .auc(A, gns_path[, j]), numeric(1))

  pos_end <- suppressWarnings(max(which(auc_path <= 0.999)))
  if (!is.finite(pos_end) || pos_end < pos_start) pos_end <- pos_start

  # 10 quantile-spaced steps from lambda.min to the undersmoothing boundary
  # (mirrors lambda2.steps <- floor(quantile(start:end, seq(0,1,.1))) in
  # Appendix S1). Collapses to the single value lambda_cv when no further
  # undersmoothing is safe (pos_end == pos_start).
  step_pos   <- unique(floor(quantile(pos_start:pos_end, probs = seq(0, 1, 0.1))))
  lambda_seq <- path_lambda[step_pos]

  best_lambda <- .collab_cv_loss(
    X, A, Y, Q1W, Q0W, folds, lambda_seq,
    ps_fit_fn = oal_ps_fn
  )

  # Cross-fitted final propensity score (bug fix #2) — out-of-fold predictions
  # (Wyss et al. Model 8) on the same shared `folds` as Q and targeting.
  g1W <- .crossfit_g(X, A, folds, best_lambda, ps_fit_fn = oal_ps_fn)

  # Cross-fitted targeting epsilon (was: .final_target(), in-sample epsilon)
  tgt <- .cv_target(Y, A, Q1W, Q0W, g1W, folds)   # CHANGED
  if (is.null(tgt))
    return(.na_est("C-TMLE1-OAL",
                   best_lambda     = best_lambda,
                   n_forced        = length(chosen_vars_out),
                   chosen_vars_out = chosen_vars_out))

  ATE_n <- mean(tgt$Q_star1 - tgt$Q_star0)
  inf <- if (cv_var && boot_var) {
    boot <- .retarget_boot_se(Y, A, Q1W, Q0W, g1W, sc$b - sc$a,
                              alpha = alpha, n_boot = n_boot)
    ATE  <- ATE_n * (sc$b - sc$a)
    list(est = ATE, se = boot$se, ci_lower = boot$ci_lower,
         ci_upper = boot$ci_upper, pvalue = .pvalue(ATE, boot$se),
         mean_IF = NA_real_, n_boot_ok = boot$n_ok)
  } else if (cv_var) {
    .ic_inference_cv(A, Y, tgt$Q_star1, tgt$Q_star0, g1W, folds,
                     ATE_n, sc$b - sc$a, alpha)
  } else {
    .ic_inference(A, Y, tgt$Q_star1, tgt$Q_star0, g1W,
                  ATE_n, sc$b - sc$a, alpha)
  }

  c(inf, list(
    estimator       = "C-TMLE1-OAL",
    best_lambda     = best_lambda,
    n_forced        = length(chosen_vars_out),
    chosen_vars_out = chosen_vars_out
  ))
}


# 9. C-TMLE1-OAL-pkg — same candidate construction as ctmle1_oal(), but the
# final selection/targeting/inference step calls the actual ctmle::ctmleGeneral()
# instead of our custom .collab_cv_loss()/.crossfit_g()/.final_target() pipeline.
#
# Purpose: a direct empirical check of ctmle1_oal() against Wyss et al.
# (2024)'s own reference implementation (their Appendix S1 calls
# ctmleGeneral(ctmletype = 1, ...) exactly this way). NOT currently wired
# into collect_results_snp() — it's a comparison tool, not (yet) a
# replacement production estimator. See compare_ctmle_oal.R for the
# side-by-side harness this was developed against.
#
# Known caveats, found by reading the installed package's source
# (CRAN ctmle v0.1.2, 2019) rather than assuming from the paper text alone:
#
#   - block_ids is accepted for interface compatibility with
#     collect_results_snp() but has no effect: ctmleGeneral() has no
#     cluster-robust SE option, only the iid varIC/varDstar formula.
#   - ctmleGeneral()'s own `alpha` argument is NOT a significance level —
#     it bounds Y away from {0,1} for the logistic-fluctuation transform.
#     The CI returned in fit$CI is hard-coded to z=1.96 internally; this
#     wrapper instead recomputes the CI from fit$var.psi using the alpha
#     this function receives, for consistency with every other estimator
#     in this file.
#   - functions_general.R::cv_general() validates a user-supplied `folds`
#     argument for length/index correctness, then unconditionally
#     overwrites it with a fresh random split a few lines later (the
#     duplicated `folds <- by(sample(1:n,n), rep(1:V, length=n), list)`
#     line). So despite passing `folds` below to align with
#     gn_candidates_cv's split, ctmleGeneral's internal best_k selection
#     does NOT actually validate against those held-out folds. Empirically
#     (see compare_ctmle_oal.R) this makes ctmleGeneral consistently prefer
#     the most-undersmoothed candidate available, while our own
#     .collab_cv_loss() — which does respect the folds it's given — more
#     often stays near lambda.min. Final ATE estimates still end up close
#     in testing so far, but this means the two methods' lambda CHOICES are
#     not validating each other; only the final targeted estimates are
#     being compared.
#
# Requires the `ctmle` package (and its dependency chain: tmle,
# SuperLearner, nnls, gam, cvAUC, ROCR, data.table) to be installed.
# Falls back to .na_est() if it isn't, rather than erroring, so sourcing
# this file doesn't require `ctmle` to be present.
 
ctmle1_oal_pkg <- function(data, alpha = 0.05, V = 10) {
 
  if (!requireNamespace("ctmle", quietly = TRUE))
    return(.na_est("C-TMLE1-OAL-pkg",
                   best_lambda     = NA_real_,
                   n_forced        = NA_integer_,
                   chosen_vars_out = character(0)))
 
  A <- data$A
  X <- .W_matrix(data)
  n <- nrow(data)
  p <- ncol(X)
 
  if (length(unique(A)) < 2)
    return(.na_est("C-TMLE1-OAL-pkg",
                   best_lambda     = NA_real_,
                   n_forced        = NA_integer_,
                   chosen_vars_out = character(0)))
 
  # Outcome LASSO / OAL penalty vector — identical to ctmle1_oal()
  qx <- .make_QX(A, X)
  pf_outcome    <- rep(1, ncol(qx$QX) - 1)
  pf_outcome[1] <- 0
 
  out_lasso <- tryCatch(
    cv.glmnet(qx$QX[, -1, drop = FALSE], data$Y, alpha = 1, nfolds = 5,
             penalty.factor = pf_outcome, standardize = TRUE),
    error = function(e) NULL
  )
  if (is.null(out_lasso))
    return(.na_est("C-TMLE1-OAL-pkg",
                   best_lambda     = NA_real_,
                   n_forced        = NA_integer_,
                   chosen_vars_out = character(0)))
 
  out_coef        <- as.numeric(coef(out_lasso, s = "lambda.min"))
  coef_names      <- rownames(coef(out_lasso, s = "lambda.min"))
  coef_w_idx      <- which(coef_names %in% colnames(X))
  chosen_vars_out <- coef_names[coef_w_idx][out_coef[coef_w_idx] != 0]
 
  pf_treatment <- rep(1, p)
  names(pf_treatment) <- colnames(X)
  pf_treatment[colnames(X) %in% chosen_vars_out] <- 0
 
  # Candidate propensity scores, with a shared fold assignment used for
  # both cv.glmnet's fit.preval AND the folds argument passed to
  # ctmleGeneral() below (see the cv_general() caveat above for why that
  # alignment is, in practice, not actually honored downstream). 
  foldid     <- caret::createFolds(A, k = V, list = FALSE)
  folds_list <- unname(split(seq_along(A), foldid))
 
  oal_fit <- tryCatch(
    cv.glmnet(X, A, family = "binomial", alpha = 1, foldid = foldid,
             penalty.factor = pf_treatment, standardize = TRUE, keep = TRUE),
    error = function(e) NULL
  )
  if (is.null(oal_fit))
    return(.na_est("C-TMLE1-OAL-pkg",
                   best_lambda     = NA_real_,
                   n_forced        = length(chosen_vars_out),
                   chosen_vars_out = chosen_vars_out))
 
  path_lambda <- oal_fit$glmnet.fit$lambda
  lambda_cv   <- oal_fit$lambda.min
  pos_start   <- which.min(abs(path_lambda - lambda_cv))
 
  gns_path    <- predict(oal_fit$glmnet.fit, newx = X, type = "response")
  gns_cv_path <- plogis(oal_fit$fit.preval[, seq_len(ncol(gns_path)), drop = FALSE])
 
  auc_path <- vapply(seq_along(path_lambda),
                     function(j) .auc(A, gns_path[, j]), numeric(1))
  pos_end  <- suppressWarnings(max(which(auc_path <= 0.999)))
  if (!is.finite(pos_end) || pos_end < pos_start) pos_end <- pos_start
 
  step_pos   <- unique(floor(quantile(pos_start:pos_end, probs = seq(0, 1, 0.1))))
  lambda_seq <- path_lambda[step_pos]
 
  gn_candidates    <- gns_path[, step_pos, drop = FALSE]
  gn_candidates_cv <- gns_cv_path[, step_pos, drop = FALSE]
 
  # ctmleGeneral() itself will do its own internal [0,1] rescaling of Y, so we pass the raw Y here.
  # Q: crude (unadjusted) arm means, per Wyss et al. Appendix S1 — they use
  # this rather than a fitted outcome model specifically so the propensity
  # model stays conservative and captures all confounder information, not
  # just what the outcome model failed to absorb. Y is passed RAW (not via
  # .scale_y()): ctmleGeneral() does its own internal [0,1] rescaling and
  # returns est/CI already back-transformed to the original Y scale.
  r0 <- rep(mean(data$Y[A == 0]), n)
  r1 <- rep(mean(data$Y[A == 1]), n)
  Q  <- cbind(r0, r1)
 
  fit <- tryCatch(
    ctmle::ctmleGeneral(
      Y = data$Y, A = A, W = as.data.frame(X), Q = Q,
      ctmletype        = 1,
      gn_candidates    = gn_candidates,
      gn_candidates_cv = gn_candidates_cv,
      family = "gaussian",
      gbound = 0.01,   # matches .clip_ps()'s [0.01, 0.99] used everywhere else
      V = V, folds = folds_list
    ),
    error = function(e) NULL
  )
  if (is.null(fit))
    return(.na_est("C-TMLE1-OAL-pkg",
                   best_lambda     = NA_real_,
                   n_forced        = length(chosen_vars_out),
                   chosen_vars_out = chosen_vars_out))
 
  est <- fit$est
  se  <- sqrt(fit$var.psi)
  z   <- qnorm(1 - alpha / 2)
 
  list(
    est             = est,
    se              = se,
    ci_lower        = est - z * se,
    ci_upper        = est + z * se,
    pvalue          = .pvalue(est, se),
    estimator       = "C-TMLE1-OAL-pkg",
    best_lambda     = lambda_seq[fit$best_k],
    n_forced        = length(chosen_vars_out),
    chosen_vars_out = chosen_vars_out
  )
}