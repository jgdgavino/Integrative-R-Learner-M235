

build_poly_basis <- function(X, degree = 2) {
  X <- as.matrix(X)
  p <- ncol(X)
  pieces <- list(X)
  if (degree >= 2) {
    for (d in 2:degree) {
      pieces[[d]] <- X^d
    }
  }
  if (p > 1 && degree >= 2) {
    idx <- combn(p, 2)
    interact <- X[, idx[1, ]] * X[, idx[2, ]]
    pieces[[length(pieces) + 1]] <- interact
  }
  B <- do.call(cbind, pieces)
  centers <- colMeans(B)
  sds <- apply(B, 2, sd)
  sds[sds == 0] <- 1
  B_std <- scale(B, center = centers, scale = sds)
  attr(B_std, "centers") <- centers
  attr(B_std, "sds") <- sds
  B_std
}

apply_poly_basis <- function(X_new, basis_ref) {
  X_new <- as.matrix(X_new)
  p <- ncol(X_new)
  pieces <- list(X_new)
  degree <- 2
  if (degree >= 2) {
    for (d in 2:degree) {
      pieces[[d]] <- X_new^d
    }
  }
  if (p > 1 && degree >= 2) {
    idx <- combn(p, 2)
    interact <- X_new[, idx[1, ]] * X_new[, idx[2, ]]
    pieces[[length(pieces) + 1]] <- interact
  }
  B <- do.call(cbind, pieces)
  centers <- attr(basis_ref, "centers")
  sds <- attr(basis_ref, "sds")
  scale(B, center = centers, scale = sds)
}

crossfit_nuisance <- function(Y, A, X, S, K = 5) {
  n <- length(Y)
  folds <- sample(rep(1:K, length.out = n))
  mu_hat <- numeric(n)
  e_hat <- numeric(n)

  for (k in 1:K) {
    train <- folds != k
    test <- folds == k

    X_aug_train <- cbind(as.matrix(X[train, ]), S[train])
    X_aug_test <- cbind(as.matrix(X[test, ]), S[test])

    fit_mu <- cv.glmnet(
      X_aug_train, Y[train],
      alpha = 0.5, nfolds = 5
    )
    mu_hat[test] <- as.numeric(
      predict(fit_mu, X_aug_test, s = "lambda.min")
    )

    fit_e <- cv.glmnet(
      X_aug_train, A[train],
      family = "binomial",
      alpha = 0.5, nfolds = 5
    )
    e_hat[test] <- as.numeric(
      predict(
        fit_e, X_aug_test,
        s = "lambda.min", type = "response"
      )
    )
  }

  e_hat <- pmax(pmin(e_hat, 0.95), 0.05)
  list(mu_hat = mu_hat, e_hat = e_hat)
}

integrative_r_loss <- function(
    theta, B_tau, B_c, Y_res, A_res, S, lambda_tau,
    lambda_c
) {
  n <- nrow(B_tau)
  p_tau <- ncol(B_tau)
  p_c <- ncol(B_c)

  alpha <- theta[1:p_tau]
  gamma <- theta[(p_tau + 1):(p_tau + p_c)]

  tau_x <- B_tau %*% alpha
  c_x <- B_c %*% gamma
  fitted <- A_res * (tau_x + (1 - S) * c_x)
  residuals <- Y_res - fitted

  loss <- mean(residuals^2) +
    lambda_tau * sum(alpha^2) +
    lambda_c * sum(gamma^2)
  loss
}

solve_integrative_rlearner <- function(
    B_tau, B_c, Y_res, A_res, S, lambda_tau, lambda_c
) {
  n <- nrow(B_tau)
  p_tau <- ncol(B_tau)
  p_c <- ncol(B_c)

  Z_tau <- A_res * B_tau
  Z_c <- A_res * (1 - S) * B_c
  Z <- cbind(Z_tau, Z_c)

  penalty <- diag(c(
    rep(lambda_tau, p_tau),
    rep(lambda_c, p_c)
  ))

  theta <- solve(
    crossprod(Z) / n + penalty,
    crossprod(Z, Y_res) / n
  )

  list(
    alpha = theta[1:p_tau],
    gamma = theta[(p_tau + 1):(p_tau + p_c)]
  )
}

cv_loss_for_lambdas <- function(
    B_tau, B_c, Y_res, A_res, S,
    lambda_tau, lambda_c, K = 5
) {
  n <- nrow(B_tau)
  folds <- sample(rep(1:K, length.out = n))
  cv_err <- numeric(K)

  for (k in 1:K) {
    train <- folds != k
    test <- folds == k

    fit <- solve_integrative_rlearner(
      B_tau[train, , drop = FALSE],
      B_c[train, , drop = FALSE],
      Y_res[train], A_res[train], S[train],
      lambda_tau, lambda_c
    )

    tau_test <- B_tau[test, , drop = FALSE] %*% fit$alpha
    c_test <- B_c[test, , drop = FALSE] %*% fit$gamma
    fitted_test <- A_res[test] *
      (tau_test + (1 - S[test]) * c_test)
    cv_err[k] <- mean((Y_res[test] - fitted_test)^2)
  }

  mean(cv_err)
}

integrative_rlearner <- function(
    data_rct, data_obs,
    outcome = "Y", treatment = "A",
    covariates = NULL,
    degree = 2,
    K_nuisance = 5,
    K_loss = 5,
    lambda_grid = NULL,
    ratio_grid = NULL,
    verbose = TRUE
) {
  if (is.null(covariates)) {
    all_cols <- colnames(data_rct)
    covariates <- setdiff(all_cols, c(outcome, treatment))
  }

  data_rct$S <- 1
  data_obs$S <- 0
  dat <- rbind(
    data_rct[, c(covariates, outcome, treatment, "S")],
    data_obs[, c(covariates, outcome, treatment, "S")]
  )

  Y <- dat[[outcome]]
  A <- dat[[treatment]]
  X <- dat[, covariates, drop = FALSE]
  S <- dat$S
  n <- nrow(dat)

  if (verbose) cat("Fitting nuisance parameters...\n")
  nuis <- crossfit_nuisance(Y, A, X, S, K = K_nuisance)
  Y_res <- Y - nuis$mu_hat
  A_res <- A - nuis$e_hat

  if (verbose) cat("Building basis functions...\n")
  B_full <- build_poly_basis(X, degree = degree)
  B_tau <- B_full
  B_c <- B_full
  p <- ncol(B_tau)

  if (is.null(lambda_grid)) {
    lambda_grid <- 10^seq(-4, 2, length.out = 15)
  }
  if (is.null(ratio_grid)) {
    ratio_grid <- c(0.01, 0.1, 0.5, 1, 2, 5, 10, 100)
  }

  if (verbose) {
    cat(
      "Grid search over", length(lambda_grid),
      "lambda_tau x", length(ratio_grid),
      "ratios...\n"
    )
  }

  best_cv <- Inf
  best_lambda_tau <- NA
  best_lambda_c <- NA

  for (lt in lambda_grid) {
    for (r in ratio_grid) {
      lc <- lt * r
      cv_err <- cv_loss_for_lambdas(
        B_tau, B_c, Y_res, A_res, S,
        lt, lc, K = K_loss
      )
      if (cv_err < best_cv) {
        best_cv <- cv_err
        best_lambda_tau <- lt
        best_lambda_c <- lc
      }
    }
  }

  if (verbose) {
    cat(
      "Optimal lambda_tau:", best_lambda_tau,
      " lambda_c:", best_lambda_c,
      " CV error:", best_cv, "\n"
    )
  }

  final_fit <- solve_integrative_rlearner(
    B_tau, B_c, Y_res, A_res, S,
    best_lambda_tau, best_lambda_c
  )

  tau_hat <- as.numeric(B_tau %*% final_fit$alpha)
  c_hat <- as.numeric(B_c %*% final_fit$gamma)

  structure(
    list(
      tau_hat = tau_hat,
      c_hat = c_hat,
      alpha = final_fit$alpha,
      gamma = final_fit$gamma,
      basis_ref = B_full,
      lambda_tau = best_lambda_tau,
      lambda_c = best_lambda_c,
      cv_error = best_cv,
      covariates = covariates,
      degree = degree,
      Y_res = Y_res,
      A_res = A_res,
      S = S,
      X = X,
      data = dat
    ),
    class = "integrative_rlearner"
  )
}

predict.integrative_rlearner <- function(object, newdata,
                                         ...) {
  X_new <- newdata[, object$covariates, drop = FALSE]
  B_new <- apply_poly_basis(X_new, object$basis_ref)
  as.numeric(B_new %*% object$alpha)
}

summary.integrative_rlearner <- function(object, ...) {
  cat("=== Integrative R-Learner Summary ===\n\n")
  cat("Sample sizes:\n")
  cat("  RCT (S=1):", sum(object$S == 1), "\n")
  cat("  OS  (S=0):", sum(object$S == 0), "\n")
  cat(
    "  Total:    ", length(object$S), "\n\n"
  )
  cat("Basis dimension:", ncol(object$basis_ref), "\n")
  cat("Optimal lambda_tau:", object$lambda_tau, "\n")
  cat("Optimal lambda_c: ", object$lambda_c, "\n")
  cat("CV error:         ", object$cv_error, "\n\n")
  cat("CATE (tau) summary:\n")
  print(summary(object$tau_hat))
  cat("\nConfounding function (c) summary:\n")
  print(summary(object$c_hat))
  invisible(object)
}
