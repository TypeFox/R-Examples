#' ZDT6 Function
#'
#' Builds and returns the two-objective ZDT6 test problem. For \eqn{m} objective it
#' is defined as follows
#' \deqn{f(\mathbf{x}) = \left(f_1(\mathbf{x}), f_2(\mathbf{x})\right)}
#' with
#' \deqn{f_1(\mathbf{x}) = 1 - \exp(-4\mathbf{x}_1)\sin^6(6\pi\mathbf{x}_1), f_2(\mathbf{x}) = g(\mathbf{x}) h(f_1(\mathbf{x}_1), g(\mathbf{x}))}
#' where
#' \deqn{g(\mathbf{x}) = 1 + 9 \left(\frac{\sum_{i = 2}^{m}\mathbf{x}_i}{m - 1}\right)^{0.25}, h(f_1, g) = 1 - \left(\frac{f_1(\mathbf{x})}{g(\mathbf{x})}\right)^2}
#' and \eqn{\mathbf{x}_i \in [0,1], i = 1, \ldots, m}.
#' This function introduced two difficulities (see reference):
#' 1. the density of solutions decreases with the closeness to the Pareto-optimal front and
#' 2. the Pareto-optimal solutions are nonuniformly distributed along the front.
#'
#' @references E. Zitzler, K. Deb, and L. Thiele. Comparison of Multiobjective
#' Evolutionary Algorithms: Empirical Results. Evolutionary Computation, 8(2):173-195, 2000
#'
#' @param dimensions [\code{integer(1)}]\cr
#'   Number of decision variables.
#' @return [\code{smoof_multi_objective_function}]
#' @export
makeZDT6Function = function(dimensions) {
  assertNumber(dimensions, lower = 2L, na.ok = FALSE)
  force(dimensions)

    # define the two-objective ZDT1 function
  fn = function(x) {
    assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
    n = length(x)
    f1 = 1 - exp(-4 * x[1]) * (sin(6 * pi * x[1]))^6
    g = 1 + 9 * (sum(x[2:n]) / (n - 1))^(0.25)
    h = 1 - (f1 / g)^2
    f2 = g * h
    return(c(f1, f2))
  }

  makeMultiObjectiveFunction(
    name = "ZDT6 Function",
    id = paste0("zdt6_", dimensions, "d_2o"),
    description = "Zitzler et al. Function N. 6",
    fn = fn,
    par.set =  makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(0, dimensions),
      upper = rep(1, dimensions),
      vector = TRUE
      ),
    n.objectives = 2L,
    ref.point = c(11, 11)
    )
}

class(makeZDT6Function) = c("function", "smoof_generator")
attr(makeZDT6Function, "name") = c("ZDT6")
attr(makeZDT6Function, "type") = c("multi-objective")
attr(makeZDT6Function, "tags") = c("multi-objective")
