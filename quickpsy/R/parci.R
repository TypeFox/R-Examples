#' Calculates the confidence intervals for the parameters
#' \code{parci} calculates the confidence intervals for the parameters
#' @keywords internal
#' @export
parci <- function(qp, ci = .95) {
  qp$parbootstrap %>% group_by_(.dots = c(qp$groups, 'parn')) %>%
    do(one_parci(., ci))
}
