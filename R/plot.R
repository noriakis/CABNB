#' Plot abundance and prevalence effects from a CAB-NB fit class object
#'
#' Plot abundance log fold changes against prevalence log odds
#' ratios and colours taxa by their FDR-significant signal type.
#'
#' @param x A `cabnb_fit` object.
#' @param alpha FDR threshold used to colour taxa.
#' @param lfc Plot `"raw"` or `"shrunken"` log fold changes.
#' @param label Label significant taxa.
#' @param label_top Maximum number of significant taxa to label.
#' @param colors Named colours for `none`, `abundance`, `prevalence`, and
#'   `abundance+prevalence`.
#' @param pch,cex Point symbol and size.
#' @param xlab,ylab,main Axis labels and title.
#' @param legend_position legend position.
#' @param ... Arguments passed to [graphics::plot.default()].
#' @return Invisibly, the plotted results with `plot_lfc` and `plot_signal`.
#' @export
plot.cabnb_fit <- function(
    x, alpha = 0.05, lfc = c("raw", "shrunken"), label = FALSE,
    label_top = 10, legend_position="bottomright",
    colors = c(none = "#777777", abundance = "#B33F62",
               prevalence = "#2F6F9F", `abundance+prevalence` = "#7A4EAB"),
    pch = 19, cex = 1.1, xlab = NULL,
    ylab = "Prevalence log odds ratio", main = NULL, ...) {
  if (!inherits(x, "cabnb_fit") || !is.data.frame(x$results)) {
    stop("x must be a cabnb_fit object with a results data frame.", call. = FALSE)
  }
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be one finite number between 0 and 1.", call. = FALSE)
  }
  lfc <- match.arg(lfc)
  lfc_col <- if (lfc == "raw") "logFC" else "shrunken_logFC"
  required <- c("taxon", lfc_col, "qvalue", "prevalence_logOR",
                "prevalence_qvalue")
  missing_cols <- setdiff(required, names(x$results))
  if (length(missing_cols)) {
    stop("cabnb_fit results are missing: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  signal_levels <- c("none", "abundance", "prevalence", "abundance+prevalence")
  if (is.null(names(colors)) || !all(signal_levels %in% names(colors))) {
    stop("colors must be named for: ", paste(signal_levels, collapse = ", "),
         call. = FALSE)
  }
  if (length(label) != 1 || is.na(label) || !is.logical(label)) {
    stop("label must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(label_top) != 1 || is.na(label_top) || label_top < 0) {
    stop("label_top must be a non-negative number.", call. = FALSE)
  }

  res <- x$results
  
  ##
  # Only the converged taxa are plotted
  ##
  
  kept <- if ("kept" %in% names(res)) !is.na(res$kept) & res$kept else
    rep(TRUE, nrow(res))
  converged <- if ("converged" %in% names(res))
    !is.na(res$converged) & res$converged else rep(TRUE, nrow(res))
  keep <- kept & is.finite(res[[lfc_col]]) &
    is.finite(res$prevalence_logOR) & converged
  d <- res[keep, , drop = FALSE]
  if (!nrow(d)) {
    stop("No fitted taxa have finite abundance and prevalence effects.", call. = FALSE)
  }
  abundance_hit <- is.finite(d$qvalue) & d$qvalue <= alpha
  prevalence_hit <- is.finite(d$prevalence_qvalue) & d$prevalence_qvalue <= alpha
  d$plot_lfc <- d[[lfc_col]]
  d$plot_signal <- ifelse(
    abundance_hit & prevalence_hit, "abundance+prevalence",
    ifelse(abundance_hit, "abundance", ifelse(prevalence_hit, "prevalence", "none"))
  )

  if (is.null(xlab)) {
    xlab <- if (lfc == "raw") "Abundance log fold change" else
      "Shrunken abundance log fold change"
  }
  graphics::plot.default(
    d$plot_lfc, d$prevalence_logOR, col = unname(colors[d$plot_signal]),
    pch = pch, cex = cex, xlab = xlab, ylab = ylab, main = main, ...
  )
  graphics::abline(h = 0, v = 0, lty = 2, col = "#999999")
  present <- signal_levels[signal_levels %in% d$plot_signal]
  legend_labels <- c(none = "Not significant", abundance = "Abundance",
                     prevalence = "Prevalence", `abundance+prevalence` = "Both")
  graphics::legend(
    legend_position, legend = legend_labels[present], col = unname(colors[present]),
    pch = pch, bty = "n", title = paste0("q <= ", format(alpha))
  )

  if (isTRUE(label)) {
    labs <- d[d$plot_signal != "none", , drop = FALSE]
    if (nrow(labs)) {
      rank_q <- pmin(labs$qvalue, labs$prevalence_qvalue, na.rm = TRUE)
      rank_q[!is.finite(rank_q)] <- Inf
      labs <- labs[order(rank_q, labs$taxon), , drop = FALSE]
      n_label <- min(nrow(labs), label_top)
      if (n_label > 0) {
        labs <- labs[seq_len(n_label), , drop = FALSE]
        graphics::text(labs$plot_lfc, labs$prevalence_logOR, labels = labs$taxon,
                       pos = 3, cex = 0.75, xpd = NA)
      }
    }
  }
  invisible(d)
}
