#' Pair comparisons of the parameters using bootstrap
#' \code{thresholdcomparisons} Calculates the bootstrap confidence intervals for the
#' difference in the parameters for two groups for all possible pairs
#' of groups
#' @keywords internal
#' @export
thresholdcomparisons <- function(qp, ci = .95) {
  one_thresholdcomparisons(qp$thresholdsbootstrap, qp$thresholds, qp$groups, ci)
}
