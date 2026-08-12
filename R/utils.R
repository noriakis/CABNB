.cabnb_check_counts <- function(counts) {
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (!is.matrix(counts) || !is.numeric(counts) || length(dim(counts)) != 2L) {
    stop("counts must be a numeric matrix with samples in rows and taxa in columns.",
         call. = FALSE)
  }
  if (!nrow(counts) || !ncol(counts)) {
    stop("counts must have at least one sample and one taxon.", call. = FALSE)
  }
  if (any(!is.finite(counts)) || any(counts < 0) ||
      any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
    stop("counts must contain finite, non-negative integer counts.", call. = FALSE)
  }
  storage.mode(counts) <- "double"
  if (is.null(rownames(counts))) rownames(counts) <- paste0("sample", seq_len(nrow(counts)))
  if (is.null(colnames(counts))) colnames(counts) <- paste0("taxon", seq_len(ncol(counts)))
  if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop("Sample and taxon names must be unique.", call. = FALSE)
  }
  counts
}

.cabnb_scalar_in_range <- function(x, name, lower, upper,
                                   lower_open = FALSE, upper_open = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (ok) {
    ok <- if (lower_open) x > lower else x >= lower
    ok <- ok && if (upper_open) x < upper else x <= upper
  }
  if (!ok) {
    left <- if (lower_open) "(" else "["
    right <- if (upper_open) ")" else "]"
    stop(name, " must be one finite number in ", left, lower, ", ", upper,
         right, ".", call. = FALSE)
  }
  invisible(x)
}

.cabnb_make_design <- function(metadata, formula, counts) {
  if (is.matrix(formula) || is.data.frame(formula)) {
    design <- as.matrix(formula)
    if (!is.numeric(design) || nrow(design) != nrow(counts)) {
      stop("A supplied design matrix must be numeric with one row per sample.",
           call. = FALSE)
    }
  } else {
    if (is.null(metadata)) {
      stop("metadata is required when formula is not a design matrix.",
           call. = FALSE)
    }
    metadata <- as.data.frame(metadata)
    if (nrow(metadata) != nrow(counts)) {
      stop("metadata must have one row per sample.", call. = FALSE)
    }
    if (!is.null(rownames(metadata)) && !identical(rownames(metadata), rownames(counts))) {
      idx <- match(rownames(counts), rownames(metadata))
      if (anyNA(idx)) {
        stop("metadata row names must match count-matrix sample names.",
             call. = FALSE)
      }
      metadata <- metadata[idx, , drop = FALSE]
    }
    design <- stats::model.matrix(formula, metadata)
  }
  if (!ncol(design) || qr(design)$rank < ncol(design)) {
    stop("The model design matrix must have full column rank.", call. = FALSE)
  }
  design
}

.cabnb_match_coef <- function(design, coef = NULL) {
  names <- colnames(design)
  if (is.null(coef)) {
    candidates <- setdiff(seq_len(ncol(design)), match("(Intercept)", names, nomatch = 0L))
    if (!length(candidates)) {
      stop("coef is required for an intercept-only model.", call. = FALSE)
    }
    return(candidates[length(candidates)])
  }
  if (length(coef) != 1L) {
    stop("coef must identify exactly one design-matrix column.", call. = FALSE)
  }
  if (is.character(coef)) {
    idx <- match(coef, names)
    if (is.na(idx)) idx <- match(make.names(coef), make.names(names))
  } else if (is.numeric(coef) && is.finite(coef) && coef == as.integer(coef)) {
    idx <- as.integer(coef)
  } else {
    idx <- NA
  }
  if (is.na(idx) || idx < 1L || idx > ncol(design)) {
    stop("coef does not identify a design-matrix column. Available columns: ",
         paste(names, collapse = ", "), call. = FALSE)
  }
  idx
}

.cabnb_geomean <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(NA)
  exp(mean(log(x)))
}

.cabnb_estimate_bias <- function(counts, library_size) {
  geometric_means <- apply(counts, 2L, .cabnb_geomean)
  usable <- is.finite(geometric_means) & geometric_means > 0
  total_factor <- rep(NA, nrow(counts))
  if (any(usable)) {
    ratios <- sweep(counts[, usable, drop = FALSE], 2L,
                    geometric_means[usable], "/")
    total_factor <- apply(ratios, 1L, function(x) {
      x <- x[is.finite(x) & x > 0]
      if (length(x)) stats::median(x) else NA
    })
  }
  depth_factor <- library_size / .cabnb_geomean(library_size)
  fallback <- depth_factor
  total_factor[!is.finite(total_factor) | total_factor <= 0] <-
    fallback[!is.finite(total_factor) | total_factor <= 0]
  total_factor <- total_factor / .cabnb_geomean(total_factor)
  log_bias <- log(total_factor / depth_factor)
  log_bias[!is.finite(log_bias)] <- 0
  log_bias <- log_bias - mean(log_bias)
  list(log_sampling_fraction = log_bias,
       sampling_fraction = exp(log_bias))
}

.cabnb_weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA)
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  x[which(cumsum(w) >= sum(w) / 2)[1L]]
}

.cabnb_shrink_dispersion <- function(alpha, mean_mu, shrink = TRUE,
                                     min_alpha = 1e-8, max_alpha = 100) {
  alpha <- pmin(pmax(alpha, min_alpha), max_alpha)
  valid <- is.finite(alpha) & is.finite(mean_mu) & mean_mu >= 0
  fallback <- stats::median(alpha[valid], na.rm = TRUE)
  if (!is.finite(fallback)) fallback <- 0.1
  alpha[!valid] <- fallback
  if (!isTRUE(shrink)) {
    return(list(estimate = alpha, trend = alpha))
  }

  log_alpha <- log(alpha)
  log_mean <- log(pmax(mean_mu, 1e-8))
  trend <- rep(stats::median(log_alpha, na.rm = TRUE), length(alpha))
  if (sum(valid) >= 8L && length(unique(round(log_mean[valid], 6))) >= 4L) {
    loess_fit <- tryCatch(
      stats::loess(log_alpha ~ log_mean, subset = valid, span = 0.75,
                   control = stats::loess.control(surface = "direct")),
      error = function(e) NULL
    )
    if (!is.null(loess_fit)) {
      predicted <- tryCatch(stats::predict(loess_fit, log_mean),
                            error = function(e) rep(NA, length(alpha)))
      trend[is.finite(predicted)] <- predicted[is.finite(predicted)]
    }
  }
  prior_var <- stats::var(log_alpha - trend, na.rm = TRUE)
  if (!is.finite(prior_var)) prior_var <- 0.10
  prior_var <- max(prior_var, 0.10)
  weight <- prior_var / (prior_var + 0.25)
  estimate <- exp(weight * log_alpha + (1 - weight) * trend)
  estimate <- pmin(pmax(estimate, min_alpha), max_alpha)
  list(estimate = estimate, trend = exp(trend))
}

.cabnb_shrink_lfc <- function(beta, se, enabled = TRUE) {
  if (!isTRUE(enabled)) {
    return(list(estimate = beta, se = se, weight = rep(1, length(beta)),
                prior_sd = NA))
  }
  ok <- is.finite(beta) & is.finite(se) & se > 0
  signal_variance <- stats::median(pmax(beta[ok]^2 - se[ok]^2, 0), na.rm = TRUE)
  if (!is.finite(signal_variance) || signal_variance <= 0) {
    signal_variance <- stats::median(beta[ok]^2, na.rm = TRUE)
  }
  if (!is.finite(signal_variance) || signal_variance <= 0) signal_variance <- 1
  weight <- rep(NA, length(beta))
  weight[ok] <- signal_variance / (signal_variance + se[ok]^2)
  list(estimate = weight * beta, se = sqrt(weight) * se,
       weight = weight, prior_sd = sqrt(signal_variance))
}

.cabnb_safe_padjust <- function(p, method = "BH") {
  out <- rep(NA, length(p))
  ok <- is.finite(p)
  if (any(ok)) out[ok] <- stats::p.adjust(p[ok], method = method)
  out
}
