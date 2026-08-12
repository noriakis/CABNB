#' Extract discoveries from a CABNB fit
#'
#' Returns either conventional abundance discoveries or discoveries from the
#' group-blind information-filtered Any test.
#'
#' @param x A `cabnb_fit` object.
#' @param type Discovery definition. `"abundance"` uses the abundance
#'   `qvalue`; `"information_filtered"` uses `information_qvalue` among
#'   information-eligible taxa.
#' @param cutoff Optional decision cutoff. By default, abundance uses `0.05`
#'   and information-filtered inference uses the calibrated cutoff stored in
#'   `x$information_filter$decision_level`.
#' @return A data frame containing the selected rows of `x$results`.
#' @examples
#' data(cabnb_example_counts)
#' data(cabnb_example_metadata)
#' fit <- cabnb_fit(
#'   cabnb_example_counts,
#'   cabnb_example_metadata,
#'   formula = ~ group,
#'   coef = "grouptreatment"
#' )
#' cabnb_discoveries(fit, type = "abundance")
#' cabnb_discoveries(fit, type = "information_filtered")
#' @export
cabnb_discoveries <- function(
    x, type = c("abundance", "information_filtered"), cutoff = NULL) {
  if (!inherits(x, "cabnb_fit") || !is.data.frame(x$results)) {
    stop("x must be a cabnb_fit object with a results data frame.",
         call. = FALSE)
  }
  type <- match.arg(type)
  results <- x$results

  if (type == "abundance") {
    if (is.null(cutoff)) cutoff <- 0.05
    .cabnb_scalar_in_range(cutoff, "cutoff", 0, 1,
                           lower_open = TRUE, upper_open = TRUE)
    selected <- results$kept & is.finite(results$qvalue) &
      results$qvalue <= cutoff
  } else {
    if (is.null(x$information_filter)) {
      stop("This fit does not contain information-filtered results.",
           call. = FALSE)
    }
    if (is.null(cutoff)) {
      cutoff <- x$information_filter$decision_level
    }
    .cabnb_scalar_in_range(cutoff, "cutoff", 0, 1,
                           lower_open = TRUE, upper_open = TRUE)
    selected <- results$information_eligible &
      is.finite(results$information_qvalue) &
      results$information_qvalue <= cutoff
  }
  results[which(selected), , drop = FALSE]
}
