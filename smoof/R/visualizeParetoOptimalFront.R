#' @title Pareto-optimal front visualization.
#'
#' @description Quickly visualize the Pareto-optimal front of a bi-criteria objective
#' function by calling the EMOA \code{\link[mco]{nsga2}} and extracting the
#' approximated Pareto-optimal front.
#'
#' @param fn [\code{smoof_multi_objective_function}]\cr
#'   Multi-objective smoof function.
#' @param ... [any]\cr
#'   Arguments passed to \code{\link[mco]{nsga2}}.
#' @examples
#' # Here we visualize the Pareto-optimal front of the bi-objective ZDT3 function
#' fn = makeZDT3Function(dimensions = 3L)
#' vis = visualizeParetoOptimalFront(fn)
#'
#' # Alternatively we can pass some more algorithm parameters to the NSGA2 algorithm
#' vis = visualizeParetoOptimalFront(fn, popsize = 1000L)
#' @return [\code{\link[ggplot2]{ggplot}}]
#' @export
visualizeParetoOptimalFront = function(fn, ...) {
  n.objectives = getNumberOfObjectives(fn)
  if (n.objectives == 1L) {
    stopf("Visualization of approximated Pareto-optimal front only possible fo multi-objective functions with two objectives at the moment.")
  }
  requirePackages("mco", why = "smoof::visualizeParetoOptimalFront")

  par.set = getParamSet(fn)

  # get approximated Pareto front
  res = mco::nsga2(fn,
    idim = getNumberOfParameters(fn),
    odim = n.objectives,
    lower.bounds = getLower(par.set),
    upper.bounds = getUpper(par.set),
    ...
  )
  eff.points = res$value[res$pareto.optimal, ]

  # transform to ggplot-friendly format
  eff.points = as.data.frame(eff.points)
  colnames(eff.points) = c("f1", "f2")

  pl = ggplot(eff.points, mapping = aes_string(x = "f1", y = "f2"))
  pl = pl + geom_point(colour = "tomato")
  pl = pl + xlab(expression(f[1])) + ylab(expression(f[2]))
  pl = pl + ggtitle(sprintf("Objective space with approximative Pareto-optimal\n
    front for the bi-criteria %s function", getName(fn)))
  return(pl)
}
