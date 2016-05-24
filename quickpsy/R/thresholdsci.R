#' Calculates the confidence intervals for the thresholds
#' \code{thresholdsci} calculates the confidence intervals for the thresholds
#' @keywords internal
#' @export
thresholdsci <- function(qp, ci = .95, method = 'percent') {
  if (length(qp$groups) == 0)
    one_thresholdsci(qp$thresholdsbootstrap, ci, method)
  else
    qp$thresholdsbootstrap %>% group_by_(.dots = qp$groups) %>%
      do(one_thresholdsci(., ci, method))
}

