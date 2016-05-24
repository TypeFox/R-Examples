#' Pair comparisons of the parameters using bootstrap
#' \code{parcomparisons} Calculates the bootstrap confidence intervals for the
#' difference in the parameters for two groups for all possible pairs
#' of groups
#' @keywords internal
#' @export

parcomparisons <- function(qp, ci = .95) {
  qp$parbootstrap %>% group_by_(.dots = 'parn') %>%
    do(one_parcomparisons(., qp$par, qp$groups, ci))
}
