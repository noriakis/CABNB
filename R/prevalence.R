#' internal function performing ridge logistic regression
#' probably replace with glmnet function
#' @noRd
.cabnb_fit_logistic_ridge <- function(y, X, lambda = 1, maxit = 50, tol = 1e-8) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  p <- ncol(X)
  beta <- rep(0, p)
  penalty <- rep(lambda, p)
  intercept_idx <- which(colnames(X) %in% c("(Intercept)", "Intercept", "X.Intercept."))
  penalty[intercept_idx] <- 0
  P <- diag(penalty, p)
  converged <- FALSE
  for (iter in seq_len(maxit)) {
    eta <- pmax(pmin(as.numeric(X %*% beta), 30), -30)
    mu <- stats::plogis(eta)
    w <- pmax(mu * (1 - mu), 1e-6)
    z <- eta + (y - mu) / w
    H <- crossprod(X, X * w) + P
    rhs <- crossprod(X, w * z)
    beta_new <- tryCatch(as.numeric(solve(H, rhs)), error = function(e) rep(NA, p))
    if (any(!is.finite(beta_new))) {
      break
    }
    if (max(abs(beta_new - beta)) < tol) {
      beta <- beta_new
      converged <- TRUE
      break
    }
    beta <- beta_new
  }
  eta <- pmax(pmin(as.numeric(X %*% beta), 30), -30)
  mu <- stats::plogis(eta)
  w <- pmax(mu * (1 - mu), 1e-6)
  H <- crossprod(X, X * w) + P
  vc <- tryCatch(solve(H), error = function(e) matrix(NA, p, p))
  se <- sqrt(diag(vc))
  list(beta = beta, se = se, converged = converged, fitted = mu)
}

.cabnb_fit_zero_taxon <- function(y, design, coef_idx, library_size,
                                  zero_depth_adjust = TRUE, glm_maxit = 50,
                                  zero_method = c("ridge", "glm"),
                                  zero_ridge_lambda = 1) {
  zero_method <- match.arg(zero_method)
  presence <- as.numeric(y > 0)
  out <- c(logOR = NA,
    SE = NA, z = NA,
    pvalue = NA, converged = NA)
  if (length(unique(presence)) < 2) {
    return(out)
  }
  zdesign <- design
  if (zero_depth_adjust) {
    log_depth <- as.numeric(scale(log(pmax(library_size, 1))))
    candidate <- cbind(zdesign, log_depth = log_depth)
    if (qr(candidate)$rank == ncol(candidate)) {
      zdesign <- candidate
    }
  }
  cname <- colnames(design)[coef_idx]
  ridx <- match(cname, colnames(zdesign))
  if (is.na(ridx)) {
    return(out)
  }
  if (zero_method == "ridge") {
    fit <- .cabnb_fit_logistic_ridge(
      presence,
      zdesign,
      lambda = zero_ridge_lambda,
      maxit = glm_maxit
    )
    beta <- fit$beta[ridx]
    se <- fit$se[ridx]
    converged <- fit$converged
  } else {
    dat <- data.frame(presence = presence, zdesign, check.names = FALSE)
    form <- stats::as.formula(paste("presence ~ 0 +", paste(colnames(zdesign), collapse = " + ")))
    fit <- tryCatch(
      suppressWarnings(stats::glm(
        form,
        data = dat,
        family = stats::binomial(),
        control = stats::glm.control(maxit = glm_maxit)
      )),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      return(out)
    }
    sm <- tryCatch(suppressWarnings(summary(fit)), error = function(e) NULL)
    if (is.null(sm) || is.null(sm$coefficients)) {
      return(out)
    }
    ridx_glm <- match(cname, rownames(sm$coefficients))
    if (is.na(ridx_glm)) {
      return(out)
    }
    beta <- sm$coefficients[ridx_glm, "Estimate"]
    se <- sm$coefficients[ridx_glm, "Std. Error"]
    converged <- isTRUE(fit$converged)
  }
  out["logOR"] <- beta
  out["SE"] <- se
  out["converged"] <- as.numeric(isTRUE(converged))
  if (!isTRUE(converged)) {
    return(out)
  }
  z <- beta / se
  out["z"] <- z
  out["pvalue"] <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  out
}

.cabnb_fit_zero_all <- function(counts, design, coef_idx, library_size,
                                zero_depth_adjust = TRUE, glm_maxit = 50,
                                zero_method = c("ridge", "glm"),
                                zero_ridge_lambda = 1) {
  zero_method <- match.arg(zero_method)
  z <- t(vapply(
    seq_len(ncol(counts)),
    function(j) .cabnb_fit_zero_taxon(
      counts[, j],
      design = design,
      coef_idx = coef_idx,
      library_size = library_size,
      zero_depth_adjust = zero_depth_adjust,
      glm_maxit = glm_maxit,
      zero_method = zero_method,
      zero_ridge_lambda = zero_ridge_lambda
    ),
    numeric(5)
  ))
  data.frame(
    taxon = colnames(counts),
    prevalence_logOR = z[, "logOR"],
    prevalence_SE = z[, "SE"],
    prevalence_z = z[, "z"],
    prevalence_pvalue = z[, "pvalue"],
    prevalence_converged = z[, "converged"] > 0,
    stringsAsFactors = FALSE
  )
}