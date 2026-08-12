#' Plot abundance and prevalence effects from a CABNB fit
#'
#' Plots abundance log fold changes against prevalence log odds ratios using base R.
#' When available, information-filtered results are used by default.
#'
#' @param x A `cabnb_fit` object.
#' @param inference Use `"information"`, conventional `"standard"` component
#'   q-values, or `"auto"` to prefer information-filtered results.
#' @param alpha Component FDR threshold used to colour taxa.
#' @param lfc Plot `"shrunken"` or `"raw"` log fold changes.
#' @param label Label discovered taxa.
#' @param label_top Maximum number of taxa to label.
#' @param colors Named colours for `none`, `abundance`, `prevalence`,
#'   `abundance+prevalence`, and `information-any`.
#' @param pch,cex Point symbol and size.
#' @param xlab,ylab,main Axis labels and title.
#' @param legend_position Legend position.
#' @param ... Arguments passed to [graphics::plot.default()].
#' @return Invisibly, the plotted result rows with `plot_lfc` and
#'   `plot_signal` columns.
#' @export
plot.cabnb_fit <- function(
    x, inference = c("auto", "information", "standard"), alpha = 0.05,
    lfc = c("shrunken", "raw"), label = FALSE, label_top = 10,
    colors = c(
      none = "#777777", abundance = "#B33F62", prevalence = "#2F6F9F",
      `abundance+prevalence` = "#7A4EAB", `information-any` = "#D89000"
    ),
    pch = 19, cex = 1.1, xlab = NULL,
    ylab = "Prevalence log odds ratio", main = NULL,
    legend_position = "bottomright", ...) {
  if (!inherits(x, "cabnb_fit") || !is.data.frame(x$results)) {
    stop("x must be a cabnb_fit object with a results data frame.",
         call. = FALSE)
  }
  .cabnb_scalar_in_range(alpha, "alpha", 0, 1,
                         lower_open = TRUE, upper_open = TRUE)
  inference <- match.arg(inference)
  if (inference == "auto") {
    inference <- if (!is.null(x$information_filter)) "information" else "standard"
  }
  if (inference == "information" && is.null(x$information_filter)) {
    stop("This fit does not contain information-filtered results.", call. = FALSE)
  }
  lfc <- match.arg(lfc)
  lfc_column <- if (lfc == "raw") "logFC" else "shrunken_logFC"
  required <- c("taxon", lfc_column, "prevalence_logOR", "converged")
  if (inference == "information") {
    required <- c(
      required, "information_reject", "information_abundance_qvalue",
      "information_prevalence_qvalue"
    )
  } else {
    required <- c(required, "qvalue", "prevalence_qvalue")
  }
  missing <- setdiff(required, names(x$results))
  if (length(missing)) {
    stop("cabnb_fit results are missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  signal_levels <- c(
    "none", "abundance", "prevalence", "abundance+prevalence",
    "information-any"
  )
  if (is.null(names(colors)) || !all(signal_levels %in% names(colors))) {
    stop("colors must be named for: ", paste(signal_levels, collapse = ", "),
         call. = FALSE)
  }
  if (!is.logical(label) || length(label) != 1L || is.na(label)) {
    stop("label must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(label_top) || length(label_top) != 1L ||
      !is.finite(label_top) || label_top < 0) {
    stop("label_top must be a non-negative number.", call. = FALSE)
  }

  results <- x$results
  plotted <- results[
    results$kept & results$converged & is.finite(results[[lfc_column]]) &
      is.finite(results$prevalence_logOR), , drop = FALSE
  ]
  if (!nrow(plotted)) {
    stop("No fitted taxa have finite abundance and prevalence effects.",
         call. = FALSE)
  }

  if (inference == "information") {
    abundance_hit <- is.finite(plotted$information_abundance_qvalue) &
      plotted$information_abundance_qvalue <= alpha
    prevalence_hit <- is.finite(plotted$information_prevalence_qvalue) &
      plotted$information_prevalence_qvalue <= alpha
    any_only <- plotted$information_reject & !abundance_hit & !prevalence_hit
    rank_q <- plotted$information_qvalue
  } else {
    abundance_hit <- is.finite(plotted$qvalue) & plotted$qvalue <= alpha
    prevalence_hit <- is.finite(plotted$prevalence_qvalue) &
      plotted$prevalence_qvalue <= alpha
    any_only <- rep(FALSE, nrow(plotted))
    rank_q <- pmin(plotted$qvalue, plotted$prevalence_qvalue, na.rm = TRUE)
  }
  plotted$plot_lfc <- plotted[[lfc_column]]
  plotted$plot_signal <- ifelse(
    abundance_hit & prevalence_hit, "abundance+prevalence",
    ifelse(abundance_hit, "abundance",
           ifelse(prevalence_hit, "prevalence",
                  ifelse(any_only, "information-any", "none")))
  )

  if (is.null(xlab)) {
    xlab <- if (lfc == "raw") "Abundance log fold change" else
      "Shrunken abundance log fold change"
  }
  if (is.null(main)) {
    main <- if (inference == "information") "CABNB information-filtered fit" else
      "CABNB fit"
  }
  graphics::plot.default(
    plotted$plot_lfc, plotted$prevalence_logOR,
    col = unname(colors[plotted$plot_signal]), pch = pch, cex = cex,
    xlab = xlab, ylab = ylab, main = main, ...
  )
  graphics::abline(h = 0, v = 0, lty = 2, col = "#999999")
  present <- signal_levels[signal_levels %in% plotted$plot_signal]
  labels <- c(
    none = "Not significant", abundance = "Abundance",
    prevalence = "Prevalence", `abundance+prevalence` = "Both",
    `information-any` = "IF"
  )
  graphics::legend(
    legend_position, legend = labels[present], col = unname(colors[present]),
    pch = pch, bty = "n"
  )

  if (isTRUE(label)) {
    labelled <- plotted[plotted$plot_signal != "none", , drop = FALSE]
    labelled$.rank_q <- rank_q[plotted$plot_signal != "none"]
    labelled$.rank_q[!is.finite(labelled$.rank_q)] <- Inf
    labelled <- labelled[order(labelled$.rank_q, labelled$taxon), , drop = FALSE]
    n_label <- min(nrow(labelled), as.integer(label_top))
    if (n_label > 0L) {
      labelled <- labelled[seq_len(n_label), , drop = FALSE]
      graphics::text(
        labelled$plot_lfc, labelled$prevalence_logOR,
        labels = labelled$taxon, pos = 3, cex = 0.75, xpd = NA
      )
    }
  }
  invisible(plotted)
}
