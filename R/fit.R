#' Fit the CABNB differential-abundance model
#'
#' Fits abundance and prevalence models to microbiome counts.
#'
#' @param counts Samples-by-taxa count matrix.
#' @param metadata Sample metadata used by `formula`.
#' @param formula Model formula; default `~ group`.
#' @param coef Coefficient name or column index to test.
#' @param library_size Sequencing depth per sample.
#' @param min_prev Minimum taxon prevalence for filtering.
#' @param min_total Minimum total taxon count for filtering.
#' @param initial_bias Bias initialization: `"median_ratio"` or `"zero"`.
#' @param maxit_comp Maximum compositional-bias iterations.
#' @param comp_tol Convergence tolerance for bias updates.
#' @param dispersion_shrink Shrink NB dispersions toward their trend.
#' @param empirical_null Apply empirical-null abundance calibration.
#' @param scale_floor Minimum abundance empirical-null scale.
#' @param alpha FDR threshold for signal classification.
#' @param zero_depth_adjust Adjust prevalence models for sequencing depth.
#' @param zero_method Prevalence method: `"ridge"` or `"glm"`.
#' @param zero_ridge_lambda Ridge penalty for prevalence models.
#' @param glm_maxit Maximum model-fitting iterations.
#' @param prevalence_empirical_null Apply empirical-null prevalence calibration.
#' @param prevalence_empirical_null_trim Tail fraction trimmed during calibration.
#' @param prevalence_empirical_scale_floor Minimum prevalence empirical-null scale.
#' @return A `cabnb_fit` object, taxon-level results are in `$results`.
#' @examples
#' data(cabnb_example_counts)
#' data(cabnb_example_metadata)
#' fit <- cabnb_fit(
#'   cabnb_example_counts,
#'   cabnb_example_metadata,
#'   formula = ~ group,
#'   coef = "grouptreatment"
#' )
#' head(fit$results)
#' @export
cabnb_fit <- function(counts, metadata, formula = ~ group,
  coef = 2, library_size = rowSums(counts), min_prev = 0.05,
  min_total = 10, initial_bias = "median_ratio",
  maxit_comp = 3, comp_tol = 1e-3, dispersion_shrink = TRUE,
  empirical_null = TRUE, scale_floor = 1.25, alpha = 0.05,
  zero_depth_adjust = TRUE,
  zero_method = c("ridge", "glm"),
  zero_ridge_lambda = 1, glm_maxit = 50,
  prevalence_empirical_null = TRUE,
  prevalence_empirical_null_trim = 0.20,
  prevalence_empirical_scale_floor = 1.25)
  {
    

  counts <- as.matrix(counts)
  design <- stats::model.matrix(formula, metadata)
  design_names <- colnames(design)
  if (length(coef) != 1) {
    stop("coef must identify exactly one design-matrix column.")
  }
  if (is.character(coef)) {
    coef_idx <- match(coef, design_names)
    if (is.na(coef_idx)) coef_idx <- match(make.names(coef), make.names(design_names))
  } else if (is.numeric(coef) && is.finite(coef) && coef == as.integer(coef)) {
    coef_idx <- as.integer(coef)
  } else {
    coef_idx <- NA
  }
  
  coef_name <- design_names[coef_idx]
  X <- design
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  comp_idx <- coef_idx

  taxa <- colnames(counts)
  
  keep <- colMeans(counts > 0) >= min_prev & colSums(counts) >= min_total
  Y <- counts[, keep, drop = FALSE]
  taxa_keep <- colnames(Y)
  
  # 1) Initial sample bias
  gm <- apply(Y, 2, function(y) exp(mean(log(y[y > 0]))))
  ok_gm <- is.finite(gm) & gm > 0
  R <- sweep(Y[, ok_gm, drop = FALSE], 2, gm[ok_gm], "/")
  total_factor <- apply(R, 1, function(x) stats::median(x[x > 0], na.rm = TRUE))
  total_factor <- total_factor / exp(mean(log(total_factor)))
  depth_factor <- library_size / exp(mean(log(library_size)))
  
  log_bias <- log(total_factor / depth_factor)
  log_bias <- log_bias - mean(log_bias)
  if (initial_bias == "zero") log_bias[] <- 0
  
  # 2) NB fitting
  fit_one <- function(y, offset, theta = NULL) {
    d <- data.frame(y = as.numeric(y), as.data.frame(X))
    d$off <- as.numeric(offset)
    
    f <- stats::as.formula(
      paste("y ~ 0 +", paste(colnames(X), collapse = " + "), "+ offset(off)")
    )
    
    fit <- if (is.null(theta)) {
      MASS::glm.nb(
        f,
        data = d,
        control = stats::glm.control(maxit = glm_maxit)
      )
    } else {
      stats::glm(
        f,
        data = d,
        family = MASS::negative.binomial(theta),
        control = stats::glm.control(maxit = glm_maxit)
      )
    }
    
    th <- if (is.null(theta)) fit$theta else theta
    
    list(
      beta = coef(fit),
      se = sqrt(diag(stats::vcov(fit))),
      mu = stats::fitted(fit),
      alpha = 1 / th,
      converged = fit$converged
    )
  }
  
  fit_all <- function(offset, theta = NULL) {
    lapply(seq_len(ncol(Y)), function(j) {
      th <- if (is.null(theta)) NULL else theta[j]
      tryCatch(
        fit_one(Y[, j], offset, th),
        error = function(e) {
          message("taxon ", j, " failed: ", e$message)
          NULL
        }
      )
    })
  }
  # 3) Iterative sample-bias update
  comp_history <- data.frame()
  
  for (it in seq_len(maxit_comp)) {
    offset <- log(library_size) + log_bias
    fits <- suppressWarnings(fit_all(offset))
    
    B <- do.call(rbind, lapply(fits, function(z) {
      if (is.null(z)) rep(NA, ncol(X)) else z$beta
    }))
    M <- sapply(fits, function(z) {
      if (is.null(z)) NA else mean(z$mu)
    })
    
    w <- sqrt(pmax(M, 1))
    ord <- order(B[, comp_idx])
    x <- B[ord, comp_idx]
    ww <- w[ord]
    ok <- is.finite(x) & is.finite(ww) & ww > 0
    
    shift <- x[ok][which(cumsum(ww[ok]) / sum(ww[ok]) >= 0.5)[1]]
    if (!is.finite(shift)) shift <- 0
    
    comp_history <- rbind(comp_history, data.frame(iteration = it, shift = shift))
    if (abs(shift) < comp_tol) break
    
    log_bias <- log_bias + as.numeric(X[, comp_idx] * shift)
    log_bias <- log_bias - mean(log_bias)
  }
  
  offset <- log(library_size) + log_bias
  
  # 4) Dispersion estimation and shrinkage
  raw <- suppressWarnings(fit_all(offset))
  raw_alpha <- sapply(raw, function(z) if (is.null(z)) NA else z$alpha)
  raw_mean <- sapply(raw, function(z) if (is.null(z)) NA else mean(z$mu))
  
  if (dispersion_shrink) {
    la <- log(raw_alpha)
    lm <- log(pmax(raw_mean, 1e-8))
    
    trend <- stats::predict(stats::loess(la ~ lm, span = 0.75),
                     newdata = data.frame(lm = lm))
    trend[!is.finite(trend)] <- stats::median(la, na.rm = TRUE)
    
    prior_var <- max(stats::var(la - trend, na.rm = TRUE), 0.05)
    sampling_var <- rep(0.25, length(la))
    w <- prior_var / (prior_var + sampling_var)
    
    disp <- exp(w * la + (1 - w) * trend)
  } else {
    disp <- raw_alpha
  }
  
  disp <- pmin(pmax(disp, 1e-8), 100)
  
  # 5) Final fixed-dispersion fit
  final <- suppressWarnings(fit_all(offset, theta = 1 / disp))
  
  beta <- sapply(final, function(z) if (is.null(z)) NA else z$beta[coef_idx])
  se <- sapply(final, function(z) if (is.null(z)) NA else z$se[coef_idx])
  z <- beta / se
  conv <- sapply(final, function(z) if (is.null(z)) NA else z$converged)
  
  # 6) Wald p-values, optionally empirical-null calibrated
  if (empirical_null) {
    emp <- .cabnb_empirical_null(
      z,
      trim = 0.2,
      scale_floor = scale_floor
    )
    p <- emp$pvalue
    center <- emp$info$center
    scale <- emp$info$scale

    # zz <- z[is.finite(z)]
    # central <- abs(zz) <= quantile(abs(zz), 0.80, na.rm = TRUE)
    
    # center <- median(zz[central], na.rm = TRUE)
    # scale <- mad(zz[central], constant = 1.4826, na.rm = TRUE)
    # scale <- max(scale, scale_floor)
    
    # z_use <- (z - center) / scale
  } else {
    center <- 0
    scale <- 1
    z_use <- z
    p <- 2 * stats::pnorm(abs(z_use), lower.tail = FALSE)

  }
  
  q <- stats::p.adjust(p, method = "BH")
  
  # 7) Prevalence model
  fit_res <- data.frame(
      taxon = taxa_keep,
      kept = TRUE,
      logFC = beta,
      log2FC = beta / log(2),
      SE = se,
      z = z,
      pvalue = p,
      qvalue = q,
      dispersion = disp,
      prevalence = colMeans(Y > 0),
      mean_count = colMeans(Y),
      signal_type = ifelse(q <= alpha, "abundance", "none"),
      converged = conv,
      stringsAsFactors = FALSE
    )
  zero_res <- .cabnb_fit_zero_all(
      Y,
      design = design,
      coef_idx = coef_idx,
      library_size = library_size,
      zero_depth_adjust = zero_depth_adjust,
      glm_maxit = glm_maxit,
      zero_method = zero_method,
      zero_ridge_lambda = zero_ridge_lambda
  )
  zero_empirical_null <- list(center = 0, scale = 1, n = 0, applied = FALSE)
  if (prevalence_empirical_null) {
      zemp <- .cabnb_empirical_null(
        zero_res$prevalence_z,
        trim = prevalence_empirical_null_trim,
        scale_floor = prevalence_empirical_scale_floor
      )
      zero_res$prevalence_analytic_pvalue <- zero_res$prevalence_pvalue
      replace <- is.finite(zemp$pvalue)
      zero_res$prevalence_pvalue[replace] <- zemp$pvalue[replace]
      zero_empirical_null <- zemp$info
  } else {
      zero_res$prevalence_analytic_pvalue <- zero_res$prevalence_pvalue
  }
  zero_res$prevalence_qvalue <- stats::p.adjust(zero_res$prevalence_pvalue,
        method = "BH")
  
  # 8) Merge the results and return
  fit_res <- merge(fit_res, zero_res, by = "taxon", all.x = TRUE, sort = FALSE)
  fit_res <- fit_res[match(taxa_keep, fit_res$taxon), , drop = FALSE]
  fit_res$prevalence_testable <- fit_res$prevalence > 0 & fit_res$prevalence < 1
  fit_res$prevalence_status <- ifelse(
    !fit_res$prevalence_testable, "no_presence_variation",
    ifelse(!is.na(fit_res$prevalence_converged) & fit_res$prevalence_converged,
           "ok", "not_converged")
  )
  abundance_hit <- is.finite(fit_res$qvalue) & fit_res$qvalue <= alpha
  prevalence_hit <- is.finite(fit_res$prevalence_qvalue) & fit_res$prevalence_qvalue <= alpha
  fit_res$signal_type <- ifelse(
    abundance_hit & prevalence_hit, "abundance+prevalence",
    ifelse(abundance_hit, "abundance", ifelse(prevalence_hit, "prevalence", "none"))
  )
  
  out <- list(
    results = fit_res,
    sample_bias = data.frame(
      sample = rownames(counts),
      library_size = library_size,
      log_sampling_fraction = log_bias,
      sampling_fraction = exp(log_bias)
    ),
    offset = offset,
    keep = keep,
    comp_history = comp_history,
    empirical_null = list(center = center, scale = scale),
    prevalence_empirical_null = zero_empirical_null,
    coefficient = coef_name,
    design = design
  )
  class(out) <- "cabnb_fit"
  return(out)
}


.cabnb_empirical_null <- function(z, trim = 0.20, scale_floor = 1, min_n = 20) {
  z <- as.numeric(z)
  ok <- is.finite(z)
  out <- rep(NA, length(z))
  info <- list(center = 0, scale = 1, n = sum(ok), applied = FALSE)
  if (sum(ok) < min_n) {
    out[ok] <- 2 * stats::pnorm(abs(z[ok]), lower.tail = FALSE)
    return(list(pvalue = out, info = info))
  }
  az <- abs(z[ok])
  cutoff <- stats::quantile(az, probs = 1 - trim, na.rm = TRUE, names = FALSE)
  central <- ok
  central[ok] <- az <= cutoff
  if (sum(central) < min_n) {
    central <- ok
  }
  center <- stats::median(z[central], na.rm = TRUE)
  scale <- stats::mad(z[central], center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    scale <- stats::sd(z[central], na.rm = TRUE)
  }
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }
  scale <- max(scale, scale_floor)
  z_cal <- (z - center) / scale
  out[ok] <- 2 * stats::pnorm(abs(z_cal[ok]), lower.tail = FALSE)
  info <- list(center = center, scale = scale, n = sum(central), applied = TRUE)
  list(pvalue = pmin(pmax(out, 0), 1), info = info)
}
