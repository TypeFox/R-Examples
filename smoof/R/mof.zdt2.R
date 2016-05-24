#' ZDT2 Function
#'
#' Builds and returns the two-objective ZDT2 test problem. The function is
#' nonconvex and resembles the ZDT1 function. For \eqn{m} objective it
#' is defined as follows
#' \deqn{f(\mathbf{x}) = \left(f_1(\mathbf{x}_1), f_2(\mathbf{x})\right)}
#' with
#' \deqn{f_1(\mathbf{x}_1) = \mathbf{x}_1, f_2(\mathbf{x}) = g(\mathbf{x}) h(f_1(\mathbf{x}_1), g(\mathbf{x}))}
#' where
#' \deqn{g(\mathbf{x}) = 1 + \frac{9}{m - 1} \sum_{i = 2}^m \mathbf{x}_i, h(f_1, g) = 1 - \left(\frac{f_1}{g}\right)^2}
#' and \eqn{\mathbf{x}_i \in [0,1], i = 1, \ldots, m}
#'
#' @references E. Zitzler, K. Deb, and L. Thiele. Comparison of Multiobjective
#' Evolutionary Algorithms: Empirical Results. Evolutionary Computation, 8(2):173-195, 2000
#'
#' @param dimensions [\code{integer(1)}]\cr
#'   Number of decision variables.
#' @return [\code{smoof_multi_objective_function}]
#' @export
makeZDT2Function = function(dimensions) {
  assertNumber(dimensions, lower = 2L, na.ok = FALSE)
  force(dimensions)

    # define the two-objective ZDT1 function
  fn = function(x) {
    assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
    n = length(x)
    f1 = x[1]
    g = 1 + 9 * sum(x[2:n]) / (n - 1)
    h = 1 - (f1 / g)^2
    f2 = g * h
    return(c(f1, f2))
  }

  makeMultiObjectiveFunction(
    name = "ZDT2 Function",
    id = paste0("zdt2_", dimensions, "d_2o"),
    description = "Zitzler et al. function N. 2",
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

class(makeZDT2Function) = c("function", "smoof_generator")
attr(makeZDT2Function, "name") = c("ZDT2")
attr(makeZDT2Function, "type") = c("multi-objective")
attr(makeZDT2Function, "tags") = c("multi-objective")
