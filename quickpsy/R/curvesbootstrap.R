#' Creates the bootstrap curves
#' \code{curvesbootstrap} creates the bootstrap curves
#' @keywords internal
#' @export
curvesbootstrap <- function(qp, xmin = NULL, xmax = NULL, log = F) {
  if (length(qp$groups) == 0)
    parboot <- qp$parbootstrap %>% group_by_('sample')
  else
    parboot <- qp$parbootstrap %>%
                group_by_(.dots = c(qp$groups, 'sample'))

  parboot %>% do(one_curve(., xmin, xmax, log, qp$groups, qp$limits,
                                   qp$psyfunguesslapses))

}



