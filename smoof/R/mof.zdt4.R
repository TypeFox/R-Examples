#' ZDT4 Function
#'
#' Builds and returns the two-objective ZDT4 test problem. For \eqn{m} objective it
#' is defined as follows
#' \deqn{f(\mathbf{x}) = \left(f_1(\mathbf{x}_1), f_2(\mathbf{x})\right)}
#' with
#' \deqn{f_1(\mathbf{x}_1) = \mathbf{x}_1, f_2(\mathbf{x}) = g(\mathbf{x}) h(f_1(\mathbf{x}_1), g(\mathbf{x}))}
#' where
#' \deqn{g(\mathbf{x}) = 1 + 10 (m - 1) + \sum_{i = 2}^{m} (\mathbf{x}_i^2 - 10\cos(4\pi\mathbf{x}_i)), h(f_1, g) = 1 - \sqrt{\frac{f_1(\mathbf{x})}{g(\mathbf{x})}}}
#' and \eqn{\mathbf{x}_i \in [0,1], i = 1, \ldots, m}.
#' This function has many Pareto-optimal fronts and is thus suited to test the
#' algorithms ability to tackle multimodal problems.
#'
#' @references E. Zitzler, K. Deb, and L. Thiele. Comparison of Multiobjective
#' Evolutionary Algorithms: Empirical Results. Evolutionary Computation, 8(2):173-195, 2000
#'
#' @param dimensions [\code{integer(1)}]\cr
#'   Number of decision variables.
#' @return [\code{smoof_multi_objective_function}]
#' @export
makeZDT4Function = function(dimensions) {
  assertNumber(dimensions, lower = 2L, na.ok = FALSE)
  force(dimensions)

  # define the two-objective ZDT1 function
  fn = function(x) {
    assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
    n = length(x)
    f1 = x[1]
    g = 1 + 10 * (n - 1) + sum(x[2:n]^2 - 10 * cos(4 * pi * x[2:n]))
    h = 1 - sqrt(f1 / g)
    f2 = g * h
    return(c(f1, f2))
  }

  makeMultiObjectiveFunction(
    name = "ZDT4 Function",
    id = paste0("zdt4_", dimensions, "d_2o"),
    description = "Zitzler et al. Function N. 4",
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

class(makeZDT4Function) = c("function", "smoof_generator")
attr(makeZDT4Function, "name") = c("ZDT4")
attr(makeZDT4Function, "type") = c("multi-objective")
attr(makeZDT4Function, "tags") = c("multi-objective")
