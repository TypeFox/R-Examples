#' Calculates bootstrap thresholds
#' \code{thresholdsbootstrap} calculates bootstrap thresholds
#' @keywords internal
#' @export
thresholdsbootstrap <- function(qp, prob = NULL, log = F) {
  if (is.null(prob)) stop('You need to specify the value of prob', call. = F)

  if (length(qp$groups) == 0)
    parboot <- qp$parbootstrap %>% group_by_('sample')
  else
    parboot <- qp$parbootstrap %>%
      group_by_(.dots = c(qp$groups, 'sample'))

  allgroups <- as.character(groups(qp$curvesbootstrap))

  parboot %>% do(one_threshold(., prob, log, allgroups,
                    qp$funname, qp$guess, qp$lapses, qp$curvesbootstrap))
}



