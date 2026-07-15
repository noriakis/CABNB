

#' Select template
.cabnb_midasim_template <- function(template) {
  template <- match.arg(template, c("throat", "ibd", "vaginal"))
  env <- new.env(parent = emptyenv())
  object_name <- switch(
    template,
    throat = "throat.otu.tab",
    ibd = "count.ibd",
    vaginal = "count.vaginal"
  )
  utils::data(list = object_name, package = "MIDASim", envir = env)
  get(object_name, envir = env, inherits = FALSE)
}

.cabnb_midasim_select_taxa <- function(counts, n_taxa) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  counts <- counts[
    rowSums(counts) > 0,
    colSums(counts > 0) > 1 & colSums(counts) > 0,
    drop = FALSE
  ]
  if (!ncol(counts)) stop("MIDASim template has no usable taxa.", call. = FALSE)
  if (n_taxa >= ncol(counts)) return(counts)
  ord <- order(colMeans(counts > 0), colMeans(counts), decreasing = TRUE)
  counts[, ord[seq_len(n_taxa)], drop = FALSE]
}

.cabnb_midasim_sample_depths <- function(template_depth, n) {
  template_depth <- template_depth[is.finite(template_depth) & template_depth > 0]
  if (!length(template_depth)) template_depth <- 5000
  pmax(100, round(sample(template_depth, n, replace = TRUE)))
}

.cabnb_pr_auc <- function(score, truth) {
  ok <- is.finite(score) & !is.na(truth)
  score <- score[ok]
  truth <- as.logical(truth[ok])
  if (!sum(truth) || !sum(!truth)) return(NA)
  curve <- PRROC::pr.curve(
    scores.class0 = score[truth],
    scores.class1 = score[!truth],
    curve = FALSE
  )
  unname(curve$auc.integral)
}

.cabnb_one_condition_metrics <- function(result, truth, target, alpha) {
  d <- merge(truth, result, by = "taxon", all.x = TRUE, sort = FALSE)
  d <- d[match(truth$taxon, d$taxon), , drop = FALSE]
  if (target == "abundance") {
    truth_positive <- d$is_abundance_da
    p_raw <- d$pvalue
    q_raw <- d$qvalue
  } else if (target == "prevalence") {
    truth_positive <- d$is_prevalence_da
    p_raw <- d$prevalence_pvalue
    q_raw <- d$prevalence_qvalue
  } else {
    truth_positive <- d$is_any_da
    p_abundance <- ifelse(is.finite(d$pvalue), d$pvalue, 1)
    p_prevalence <- ifelse(is.finite(d$prevalence_pvalue), d$prevalence_pvalue, 1)
    p_raw <- pmin(1, 2 * pmin(p_abundance, p_prevalence))
    q_raw <- stats::p.adjust(p_raw, method = "BH")
  }
  truth_positive <- as.logical(truth_positive)
  tested <- is.finite(p_raw)
  p_eval <- ifelse(tested, p_raw, 1)
  q_eval <- ifelse(is.finite(q_raw), q_raw, 1)
  called <- q_eval <= alpha
  tp <- sum(called & truth_positive)
  fp <- sum(called & !truth_positive)
  fn <- sum(!called & truth_positive)
  tn <- sum(!called & !truth_positive)
  score <- -log10(pmax(p_eval, .Machine$double.xmin))
  data.frame(
    target = target,
    n_taxa = nrow(d),
    n_truth = sum(truth_positive),
    n_tested = sum(tested),
    discoveries = sum(called),
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn,
    fdr = fp / max(tp + fp, 1),
    power = tp / max(tp + fn, 1),
    recall = tp / max(tp + fn, 1),
    precision = tp / max(tp + fp, 1),
    specificity = tn / max(tn + fp, 1),
    pr_auc = .cabnb_pr_auc(score, truth_positive),
    stringsAsFactors = FALSE
  )
}

#' Run and evaluate CABNB for one MIDASim condition
#'
#' @param template MIDASim template: `"throat"`, `"ibd"`, or `"vaginal"`.
#' @param n_per_group Number of control and treatment samples per group.
#' @param n_taxa Number of template taxa retained for simulation.
#' @param da_prop Proportion of taxa receiving abundance effects.
#' @param prevalence_prop Proportion receiving prevalence effects.
#' @param effect_log_fc Absolute natural-log abundance effect size.
#' @param prevalence_drop Probability that a non-zero treatment count is set to
#'   zero for a prevalence-effect taxon.
#' @param alpha FDR threshold for discoveries.
#' @param seed Random seed for this condition.
#' @param cabnb_args Named list of additional arguments passed to cabnb_fit().
#' @param save_prefix Optional file prefix. When supplied, CSV files for metrics,
#'   taxon results, and truth plus one RDS file are written.
#' @return A list containing settings, counts, metadata, truth, fit, and metrics.
cabnb_run_midasim_one_condition <- function(
    template = c("throat", "ibd", "vaginal"), n_per_group = 40, n_taxa = 80,
    da_prop = 0.10, prevalence_prop = 0.05, effect_log_fc = log(2), prevalence_drop = 0.35,
    alpha = 0.05, seed = 1, cabnb_args = list(), save_prefix = NULL) {
  if (!requireNamespace("MIDASim", quietly = TRUE)) {
    stop("MIDASim is required. Install it before running this script.", call. = FALSE)
  }
  if (!requireNamespace("PRROC", quietly = TRUE)) {
    stop("PRROC is required for PR-AUC calculation.", call. = FALSE)
  }
  if (!exists("cabnb_fit", mode = "function")) {
    stop("Load CABNB before calling cabnb_run_midasim_one_condition().", call. = FALSE)
  }
  template <- match.arg(template)
  stopifnot(n_per_group >= 2, n_taxa >= 2, da_prop > 0, da_prop < 1,
            prevalence_prop >= 0, prevalence_prop < 1, alpha > 0, alpha < 1)
  set.seed(seed)

  template_counts <- .cabnb_midasim_select_taxa(
    .cabnb_midasim_template(template), n_taxa
  )
  template_depth <- rowSums(template_counts)
  setup <- MIDASim::MIDASim.setup(template_counts)
  n_taxa_actual <- setup$n.taxa
  taxa <- setup$taxa.names
  n <- 2 * n_per_group
  group <- factor(rep(c("control", "treatment"), each = n_per_group),
                  levels = c("control", "treatment"))
  group01 <- as.numeric(group == "treatment")

  n_da <- min(n_taxa_actual, max(1, round(da_prop * n_taxa_actual)))
  n_prev <- min(n_taxa_actual - n_da,
                max(0, round(prevalence_prop * n_taxa_actual)))
  abundance_idx <- sample(seq_len(n_taxa_actual), n_da)
  remaining <- setdiff(seq_len(n_taxa_actual), abundance_idx)
  prevalence_idx <- if (n_prev) sample(remaining, n_prev) else integer()
  true_log_fc <- numeric(n_taxa_actual)
  true_log_fc[abundance_idx] <- sample(rep(c(-1, 1), length.out = n_da)) *
    effect_log_fc

  base_rel <- setup$mean.rel.abund / sum(setup$mean.rel.abund)
  treatment_rel <- base_rel * exp(true_log_fc)
  treatment_rel <- treatment_rel / sum(treatment_rel)
  planned_depth <- .cabnb_midasim_sample_depths(template_depth, n)

  taxa_presence <- setup$taxa.1.prop
  sample_presence <- rep(mean(taxa_presence), n_per_group)
  control_model <- MIDASim::MIDASim.modify(
    setup,
    lib.size = planned_depth[seq_len(n_per_group)],
    mean.rel.abund = base_rel,
    sample.1.prop = sample_presence,
    taxa.1.prop = taxa_presence
  )
  treatment_model <- MIDASim::MIDASim.modify(
    setup,
    lib.size = planned_depth[n_per_group + seq_len(n_per_group)],
    mean.rel.abund = treatment_rel,
    sample.1.prop = sample_presence,
    taxa.1.prop = taxa_presence
  )
  control_counts <- MIDASim::MIDASim(control_model)$sim_count
  treatment_counts <- MIDASim::MIDASim(treatment_model)$sim_count
  counts <- rbind(control_counts, treatment_counts)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  rownames(counts) <- sprintf("Sample_%03d", seq_len(n))
  colnames(counts) <- taxa

  prevalence_effect <- logical(n_taxa_actual)
  treatment_rows <- which(group01 == 1)
  for (j in prevalence_idx) {
    nonzero <- treatment_rows[counts[treatment_rows, j] > 0]
    if (length(nonzero)) {
      counts[nonzero[stats::runif(length(nonzero)) < prevalence_drop], j] <- 0
    }
    prevalence_effect[j] <- TRUE
  }

  metadata <- data.frame(
    group = group,
    planned_library_size = planned_depth,
    observed_library_size = rowSums(counts),
    row.names = rownames(counts)
  )
  formula <- ~ group
  truth <- data.frame(
    taxon = taxa,
    is_abundance_da = true_log_fc != 0,
    is_prevalence_da = prevalence_effect,
    true_log_fc = true_log_fc,
    stringsAsFactors = FALSE
  )
  truth$is_any_da <- truth$is_abundance_da | truth$is_prevalence_da

  fit_call <- utils::modifyList(
    list(counts = counts, metadata = metadata, formula = formula,
         coef = "grouptreatment", alpha = alpha),
    cabnb_args
  )
  fit <- do.call(cabnb_fit, fit_call)
  metrics <- do.call(rbind, lapply(
    c("abundance", "prevalence", "any"),
    function(target) .cabnb_one_condition_metrics(
      fit$results, truth, target = target, alpha = alpha
    )
  ))
  rownames(metrics) <- NULL
  settings <- list(
    template = template, n_per_group = n_per_group, n_taxa = n_taxa_actual,
    da_prop = da_prop, prevalence_prop = prevalence_prop,
    effect_log_fc = effect_log_fc, prevalence_drop = prevalence_drop,
    alpha = alpha, seed = seed
  )
  out <- list(settings = settings, counts = counts, metadata = metadata,
              truth = truth, fit = fit, metrics = metrics)

  if (!is.null(save_prefix)) {
    parent <- dirname(save_prefix)
    if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
    utils::write.csv(metrics, paste0(save_prefix, "_metrics.csv"), row.names = FALSE)
    utils::write.csv(fit$results, paste0(save_prefix, "_taxon_results.csv"),
                     row.names = FALSE)
    utils::write.csv(truth, paste0(save_prefix, "_truth.csv"), row.names = FALSE)
    saveRDS(out, paste0(save_prefix, ".rds"))
  }
  out
}
