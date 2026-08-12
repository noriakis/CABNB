.cabnb_any_min_p <- function(abundance_p, prevalence_p) {
  components <- cbind(as.numeric(abundance_p), as.numeric(prevalence_p))
  n_components <- rowSums(is.finite(components))
  min_p <- apply(components, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) min(x) else NA
  })
  out <- rep(NA, nrow(components))
  testable <- n_components > 0L
  out[testable] <- pmin(1, n_components[testable] * min_p[testable])
  out
}

.cabnb_apply_information_filter <- function(
    fit, counts, min_effective_n = 10, max_sample_fraction = 0.35,
    min_prevalence = 0, decision_level = 0.20,
    p_adjust_method = "BH") {
  .cabnb_scalar_in_range(min_effective_n, "information_min_effective_n", 1, Inf)
  .cabnb_scalar_in_range(max_sample_fraction,
                         "information_max_sample_fraction", 0, 1,
                         lower_open = TRUE)
  .cabnb_scalar_in_range(min_prevalence, "information_min_prevalence", 0, 1)
  .cabnb_scalar_in_range(decision_level, "information_decision_level", 0, 1,
                         lower_open = TRUE)

  results <- fit$results
  idx <- match(results$taxon, colnames(counts))
  total <- colSums(counts)[idx]
  sum_squares <- colSums(counts^2)[idx]
  maximum <- apply(counts, 2L, max)[idx]
  pooled_prevalence <- colMeans(counts > 0)[idx]
  effective_n <- total^2 / pmax(sum_squares, 1)
  max_fraction <- maximum / pmax(total, 1)

  any_pvalue <- .cabnb_any_min_p(
    results$pvalue, results$prevalence_pvalue
  )
  information_ok <- results$kept & is.finite(effective_n) &
    is.finite(max_fraction) & is.finite(pooled_prevalence) &
    effective_n >= min_effective_n &
    max_fraction <= max_sample_fraction &
    pooled_prevalence >= min_prevalence
  any_eligible <- information_ok & is.finite(any_pvalue)
  abundance_eligible <- information_ok & is.finite(results$pvalue)
  prevalence_eligible <- information_ok &
    is.finite(results$prevalence_pvalue)

  any_qvalue <- rep(NA, nrow(results))
  any_qvalue[any_eligible] <- .cabnb_safe_padjust(
    any_pvalue[any_eligible], p_adjust_method
  )
  abundance_qvalue <- rep(NA, nrow(results))
  abundance_qvalue[abundance_eligible] <- .cabnb_safe_padjust(
    results$pvalue[abundance_eligible], p_adjust_method
  )
  prevalence_qvalue <- rep(NA, nrow(results))
  prevalence_qvalue[prevalence_eligible] <- .cabnb_safe_padjust(
    results$prevalence_pvalue[prevalence_eligible], p_adjust_method
  )

  results$any_pvalue <- any_pvalue
  results$any_qvalue_unfiltered <- .cabnb_safe_padjust(
    any_pvalue, p_adjust_method
  )
  results$information_effective_n <- effective_n
  results$information_max_sample_fraction <- max_fraction
  results$information_pooled_prevalence <- pooled_prevalence
  results$information_eligible <- any_eligible
  results$information_qvalue <- any_qvalue
  results$information_reject <- is.finite(any_qvalue) &
    any_qvalue <= decision_level
  results$information_abundance_eligible <- abundance_eligible
  results$information_abundance_qvalue <- abundance_qvalue
  results$information_abundance_reject <- is.finite(abundance_qvalue) &
    abundance_qvalue <= 0.05
  results$information_prevalence_eligible <- prevalence_eligible
  results$information_prevalence_qvalue <- prevalence_qvalue
  results$information_prevalence_reject <- is.finite(prevalence_qvalue) &
    prevalence_qvalue <= 0.05

  # Short aliases retained for downstream result-table workflows.
  results$p_any <- any_pvalue
  results$q_any <- any_qvalue
  fit$results <- results
  fit$information_filter <- list(
    group_label_independent = TRUE,
    min_effective_n = min_effective_n,
    max_sample_fraction = max_sample_fraction,
    min_prevalence = min_prevalence,
    p_adjust_method = p_adjust_method,
    decision_level = decision_level,
    decision_level_is_conventional_fdr = FALSE,
    n_eligible = sum(any_eligible),
    n_rejected = sum(results$information_reject),
    abundance_decision_level = 0.05,
    abundance_decision_level_is_conventional_fdr = TRUE,
    n_abundance_eligible = sum(abundance_eligible),
    n_abundance_rejected = sum(results$information_abundance_reject),
    n_prevalence_eligible = sum(prevalence_eligible),
    n_prevalence_rejected = sum(results$information_prevalence_reject)
  )
  fit
}
