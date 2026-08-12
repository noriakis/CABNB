#' Fit the CABNB model
#'
#' Fits a compositionality-aware negative-binomial abundance model, an
#' optional prevalence model, and a group-blind information filter. Counts
#' must have samples in rows and taxa in columns.
#'
#' @param counts Numeric samples-by-taxa count matrix.
#' @param metadata Sample metadata used by `formula`.
#' @param formula Model formula; default `~ group`.
#' @param coef Coefficient name or design-matrix column index to test. By
#'   default, the last non-intercept coefficient is used.
#' @param library_size Optional sequencing depth per sample.
#'   Defaults to row sums of `counts`.
#' @param min_prevalence Minimum taxon prevalence for model fitting.
#' @param min_total_count Minimum total taxon count for model fitting.
#' @param min_mean_count Minimum mean taxon count for model fitting.
#' @param initial_bias Bias initialization: `"median_ratio"` or `"zero"`.
#' @param maxit_comp Maximum compositional-bias iterations.
#' @param comp_tol Convergence tolerance for compositional-bias updates.
#' @param dispersion_shrink Shrink NB dispersions toward their mean trend.
#' @param lfc_shrink Apply normal-prior shrinkage to reported log fold changes.
#' @param empirical_null Apply robust empirical-null calibration to abundance
#'   Wald statistics.
#' @param empirical_null_trim Upper tail fraction excluded from calibration.
#' @param empirical_null_scale_floor Minimum abundance empirical-null scale.
#' @param p_adjust_method Multiple-testing adjustment method.
#' @param alpha FDR threshold used for component signal classification.
#' @param zero_model Fit the prevalence model on presence/absence indicators.
#' @param zero_depth_adjust Adjust prevalence models for sequencing depth.
#' @param zero_method Prevalence method: `"ridge"` or `"glm"`.
#' @param zero_ridge_lambda Ridge penalty for prevalence models.
#' @param glm_maxit Maximum model-fitting iterations.
#' @param prevalence_empirical_null Apply robust empirical-null prevalence
#'   calibration.
#' @param prevalence_empirical_scale_floor Minimum prevalence empirical-null
#'   scale.
#' @param information_filter Apply the group-blind information filter.
#' @param information_min_effective_n Minimum pooled count-support effective
#'   sample size, `(sum(y)^2 / sum(y^2))`.
#' @param information_max_sample_fraction Maximum fraction of a taxon's total
#'   count contributed by one sample.
#' @param information_min_prevalence Minimum pooled prevalence.
#' @param information_decision_level Calibrated cutoff for the filtered Any
#'   q-value. This is not a conventional FDR guarantee.
#' @param verbose Report taxa that fail model fitting.
#' @return A `cabnb_fit` object. Taxon-level results are in `$results`.
#' @examples
#' data(cabnb_example_counts)
#' data(cabnb_example_metadata)
#' fit <- cabnb_fit(
#'   cabnb_example_counts,
#'   cabnb_example_metadata,
#'   formula = ~ group,
#'   coef = "grouptreatment"
#' )
#' fit
#' head(fit$results)
#' @export
cabnb_fit <- function(
    counts, metadata = NULL, formula = ~ group, coef = NULL,
    library_size = NULL, min_prevalence = 0.05, min_total_count = 10,
    min_mean_count = 0, initial_bias = c("median_ratio", "zero"),
    maxit_comp = 3, comp_tol = 1e-3, dispersion_shrink = FALSE,
    lfc_shrink = TRUE, empirical_null = TRUE, empirical_null_trim = 0.20,
    empirical_null_scale_floor = 1.25, p_adjust_method = "BH", alpha = 0.05,
    zero_model = TRUE, zero_depth_adjust = TRUE,
    zero_method = c("ridge", "glm"), zero_ridge_lambda = 1,
    glm_maxit = 50, prevalence_empirical_null = TRUE,
    prevalence_empirical_scale_floor = 1.25,
    information_filter = TRUE, information_min_effective_n = 10,
    information_max_sample_fraction = 0.35, information_min_prevalence = 0,
    information_decision_level = 0.20, verbose = FALSE) {
  counts <- .cabnb_check_counts(counts)
  initial_bias <- match.arg(initial_bias)
  zero_method <- match.arg(zero_method)

  .cabnb_scalar_in_range(alpha, "alpha", lower = 0, upper = 1,
                         lower_open = TRUE, upper_open = TRUE)
  .cabnb_scalar_in_range(min_prevalence, "min_prevalence", 0, 1)
  .cabnb_scalar_in_range(min_total_count, "min_total_count", 0, Inf)
  .cabnb_scalar_in_range(min_mean_count, "min_mean_count", 0, Inf)
  .cabnb_scalar_in_range(maxit_comp, "maxit_comp", 0, Inf)
  .cabnb_scalar_in_range(comp_tol, "comp_tol", 0, Inf)

  if (is.null(library_size)) {
    library_size <- rowSums(counts)
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != nrow(counts) ||
      any(!is.finite(library_size) | library_size <= 0)) {
    stop("library_size must contain one positive finite value per sample.",
         call. = FALSE)
  }

  design <- .cabnb_make_design(metadata, formula, counts)
  coef_idx <- .cabnb_match_coef(design, coef)
  coef_name <- colnames(design)[coef_idx]
  fit_design <- design
  colnames(fit_design) <- paste0("x", seq_len(ncol(fit_design)))

  keep <- colMeans(counts > 0) >= min_prevalence &
    colSums(counts) >= min_total_count &
    colMeans(counts) >= min_mean_count
  if (!any(keep)) {
    stop("No taxa passed filtering. Relax the minimum count/prevalence settings.",
         call. = FALSE)
  }
  y_fit <- counts[, keep, drop = FALSE]

  bias <- .cabnb_estimate_bias(y_fit, library_size)
  median_ratio_log_bias <- bias$log_sampling_fraction
  log_bias <- if (initial_bias == "zero") {
    rep(0, nrow(counts))
  } else {
    median_ratio_log_bias
  }

  fit_one <- function(y, offset, theta = NULL) {
    dat <- data.frame(y = as.numeric(y), as.data.frame(fit_design))
    dat$.offset <- as.numeric(offset)
    form <- stats::as.formula(paste(
      "y ~ 0 +", paste(colnames(fit_design), collapse = " + "),
      "+ offset(.offset)"
    ))
    fit <- if (is.null(theta)) {
      MASS::glm.nb(form, data = dat,
                   control = stats::glm.control(maxit = glm_maxit))
    } else {
      stats::glm(form, data = dat, family = MASS::negative.binomial(theta),
                 control = stats::glm.control(maxit = glm_maxit))
    }
    theta_fit <- if (is.null(theta)) fit$theta else theta
    list(
      beta = stats::coef(fit),
      se = sqrt(diag(stats::vcov(fit))),
      mu = stats::fitted(fit),
      dispersion = 1 / theta_fit,
      converged = isTRUE(fit$converged),
      message = NULL
    )
  }

  fit_all <- function(offset, theta = NULL) {
    fits <- lapply(seq_len(ncol(y_fit)), function(j) {
      theta_j <- if (is.null(theta)) NULL else theta[j]
      tryCatch(
        suppressWarnings(fit_one(y_fit[, j], offset, theta_j)),
        error = function(e) list(
          beta = rep(NA, ncol(fit_design)),
          se = rep(NA, ncol(fit_design)),
          mu = rep(NA, nrow(y_fit)),
          dispersion = NA, converged = FALSE,
          message = conditionMessage(e)
        )
      )
    })
    names(fits) <- colnames(y_fit)
    if (isTRUE(verbose)) {
      failed <- vapply(fits, function(x) !isTRUE(x$converged), logical(1))
      if (any(failed)) {
        message("CABNB: ", sum(failed), " taxon model(s) did not converge.")
      }
    }
    fits
  }

  comp_history <- data.frame(
    iteration = integer(), shift = numeric(), stringsAsFactors = FALSE
  )
  if (maxit_comp > 0) {
    for (iteration in seq_len(as.integer(maxit_comp))) {
      offset_iteration <- log(library_size) + log_bias
      fits_iteration <- fit_all(offset_iteration)
      beta_iteration <- vapply(
        fits_iteration, function(x) x$beta[coef_idx], numeric(1)
      )
      mean_iteration <- vapply(
        fits_iteration, function(x) mean(x$mu, na.rm = TRUE), numeric(1)
      )
      shift <- .cabnb_weighted_median(
        beta_iteration, sqrt(pmax(mean_iteration, 1))
      )
      if (!is.finite(shift)) shift <- 0
      comp_history <- rbind(
        comp_history,
        data.frame(iteration = iteration, shift = shift)
      )
      if (abs(shift) < comp_tol) break
      log_bias <- log_bias + as.numeric(fit_design[, coef_idx] * shift)
      log_bias <- log_bias - mean(log_bias)
    }
  }

  offset <- log(library_size) + log_bias
  raw_fits <- fit_all(offset)
  raw_dispersion <- vapply(raw_fits, function(x) x$dispersion, numeric(1))
  raw_mean <- vapply(raw_fits, function(x) mean(x$mu, na.rm = TRUE), numeric(1))
  dispersion <- .cabnb_shrink_dispersion(
    raw_dispersion, raw_mean, shrink = dispersion_shrink
  )
  final_fits <- fit_all(offset, theta = 1 / dispersion$estimate)

  beta <- vapply(final_fits, function(x) x$beta[coef_idx], numeric(1))
  se <- vapply(final_fits, function(x) x$se[coef_idx], numeric(1))
  converged <- vapply(final_fits, function(x) isTRUE(x$converged), logical(1))
  z <- beta / se
  z[!converged | !is.finite(z)] <- NA

  if (isTRUE(empirical_null)) {
    abundance_null <- .cabnb_empirical_null(
      z, trim = empirical_null_trim,
      scale_floor = empirical_null_scale_floor
    )
    pvalue <- abundance_null$pvalue
  } else {
    abundance_null <- list(
      info = list(center = 0, scale = 1, n = sum(is.finite(z)), applied = FALSE)
    )
    pvalue <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  }
  qvalue <- .cabnb_safe_padjust(pvalue, p_adjust_method)
  shrink <- .cabnb_shrink_lfc(beta, se, enabled = lfc_shrink)

  fit_results <- data.frame(
    taxon = colnames(y_fit), kept = TRUE, filtered = FALSE,
    logFC = beta, log2FC = beta / log(2), SE = se, z = z,
    pvalue = pvalue, qvalue = qvalue,
    shrunken_logFC = shrink$estimate,
    shrunken_log2FC = shrink$estimate / log(2),
    shrunken_SE = shrink$se, lfc_shrink_weight = shrink$weight,
    dispersion = dispersion$estimate,
    dispersion_raw = raw_dispersion,
    dispersion_trend = dispersion$trend,
    prevalence = colMeans(y_fit > 0), mean_count = colMeans(y_fit),
    converged = converged, stringsAsFactors = FALSE
  )

  if (isTRUE(zero_model)) {
    prevalence_results <- .cabnb_fit_zero_all(
      y_fit, design = design, coef_idx = coef_idx,
      library_size = library_size, zero_depth_adjust = zero_depth_adjust,
      glm_maxit = glm_maxit, zero_method = zero_method,
      zero_ridge_lambda = zero_ridge_lambda
    )
    prevalence_null <- list(center = 0, scale = 1, n = 0, applied = FALSE)
    prevalence_results$prevalence_analytic_pvalue <-
      prevalence_results$prevalence_pvalue
    if (isTRUE(prevalence_empirical_null)) {
      calibrated <- .cabnb_empirical_null(
        prevalence_results$prevalence_z,
        trim = empirical_null_trim,
        scale_floor = prevalence_empirical_scale_floor
      )
      replace <- is.finite(calibrated$pvalue)
      prevalence_results$prevalence_pvalue[replace] <-
        calibrated$pvalue[replace]
      prevalence_null <- calibrated$info
    }
    prevalence_results$prevalence_qvalue <- .cabnb_safe_padjust(
      prevalence_results$prevalence_pvalue, p_adjust_method
    )
    fit_results <- merge(
      fit_results, prevalence_results, by = "taxon", all.x = TRUE, sort = FALSE
    )
    fit_results <- fit_results[
      match(colnames(y_fit), fit_results$taxon), , drop = FALSE
    ]
  } else {
    prevalence_null <- list(center = 0, scale = 1, n = 0, applied = FALSE)
    fit_results$prevalence_logOR <- NA
    fit_results$prevalence_SE <- NA
    fit_results$prevalence_z <- NA
    fit_results$prevalence_analytic_pvalue <- NA
    fit_results$prevalence_pvalue <- NA
    fit_results$prevalence_qvalue <- NA
    fit_results$prevalence_converged <- NA
  }

  abundance_hit <- is.finite(fit_results$qvalue) &
    fit_results$qvalue <= alpha
  prevalence_hit <- is.finite(fit_results$prevalence_qvalue) &
    fit_results$prevalence_qvalue <= alpha
  fit_results$signal_type <- ifelse(
    abundance_hit & prevalence_hit, "abundance+prevalence",
    ifelse(abundance_hit, "abundance",
           ifelse(prevalence_hit, "prevalence", "none"))
  )

  results <- data.frame(
    taxon = colnames(counts), kept = keep, filtered = !keep,
    stringsAsFactors = FALSE
  )
  results <- merge(results, fit_results, by = c("taxon", "kept", "filtered"),
                   all.x = TRUE, sort = FALSE)
  results <- results[match(colnames(counts), results$taxon), , drop = FALSE]

  fit <- list(
    results = results,
    coefficient = coef_name,
    design = design,
    metadata = metadata,
    sample_bias = data.frame(
      sample = rownames(counts), library_size = library_size,
      median_ratio_log_sampling_fraction = median_ratio_log_bias,
      log_sampling_fraction = log_bias,
      sampling_fraction = exp(log_bias), stringsAsFactors = FALSE
    ),
    offset = offset, keep = keep, comp_history = comp_history,
    empirical_null = abundance_null$info,
    prevalence_empirical_null = prevalence_null,
    lfc_prior_sd = shrink$prior_sd,
    initial_bias = initial_bias,
    alpha = alpha,
    p_adjust_method = p_adjust_method,
    call = match.call()
  )

  if (isTRUE(information_filter)) {
    fit <- .cabnb_apply_information_filter(
      fit, counts = counts,
      min_effective_n = information_min_effective_n,
      max_sample_fraction = information_max_sample_fraction,
      min_prevalence = information_min_prevalence,
      decision_level = information_decision_level,
      p_adjust_method = p_adjust_method
    )
  } else {
    fit$information_filter <- NULL
  }

  class(fit) <- "cabnb_fit"
  return(fit)
}

#' @export
print.cabnb_fit <- function(x, ...) {
  results <- x$results
  alpha <- if (!is.null(x$alpha)) x$alpha else 0.05
  cat("CABNB fit\n")
  cat("  coefficient:", x$coefficient, "\n")
  cat("  fitted taxa:", sum(results$kept), "of", nrow(results), "\n")
  cat("  abundance discoveries at q <=", alpha, ":",
      sum(is.finite(results$qvalue) & results$qvalue <= alpha), "\n")
  cat("  prevalence discoveries at q <=", alpha, ":",
      sum(is.finite(results$prevalence_qvalue) &
            results$prevalence_qvalue <= alpha), "\n")
  if (!is.null(x$information_filter)) {
    info <- x$information_filter
    cat("  information-eligible taxa:", info$n_eligible, "\n")
    cat("  information-filtered Any discoveries at calibrated cutoff",
        info$decision_level, ":", info$n_rejected, "\n")
  }
  invisible(x)
}

.cabnb_empirical_null <- function(z, trim = 0.20, scale_floor = 1, min_n = 20) {
  z <- as.numeric(z)
  ok <- is.finite(z)
  pvalue <- rep(NA, length(z))
  info <- list(center = 0, scale = 1, n = sum(ok), applied = FALSE)
  if (sum(ok) < min_n) {
    pvalue[ok] <- 2 * stats::pnorm(abs(z[ok]), lower.tail = FALSE)
    return(list(pvalue = pvalue, info = info))
  }
  cutoff <- stats::quantile(
    abs(z[ok]), probs = 1 - trim, na.rm = TRUE, names = FALSE
  )
  central <- ok
  central[ok] <- abs(z[ok]) <= cutoff
  if (sum(central) < min_n) central <- ok
  center <- stats::median(z[central], na.rm = TRUE)
  scale <- stats::mad(
    z[central], center = center, constant = 1.4826, na.rm = TRUE
  )
  if (!is.finite(scale) || scale <= 0) {
    scale <- stats::sd(z[central], na.rm = TRUE)
  }
  if (!is.finite(scale) || scale <= 0) scale <- 1
  scale <- max(scale, scale_floor)
  pvalue[ok] <- 2 * stats::pnorm(
    abs((z[ok] - center) / scale), lower.tail = FALSE
  )
  list(
    pvalue = pmin(pmax(pvalue, 0), 1),
    info = list(center = center, scale = scale, n = sum(central), applied = TRUE)
  )
}
