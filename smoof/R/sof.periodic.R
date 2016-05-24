#' Periodic Function
#'
#' Single-objective two-dimensional test function. The formula is given as
#' \deqn{f(\mathbf{x}) = 1 + \sin^2(\mathbf{x}_1) + \sin^2(\mathbf{x}_2) - 0.1e^{-(\mathbf{x}_1^2 + \mathbf{x}_2^2)}}
#' subject to the constraints \eqn{\mathbf{x}_i \in [-10, 10], i = 1, 2.}
#'
#' @references M. M. Ali, C. Khompatraporn, Z. B. Zabinsky, A Numerical
#' Evaluation of Several Stochastic Algorithms on Selected Continuous Global
#' Optimization Test Problems, Journal of Global Optimization, vol. 31, pp.
#' 635-672, 2005.
#'
#' @template ret_smoof_single
#' @export
makePeriodicFunction = function() {
  makeSingleObjectiveFunction(
    name = "Periodic Function",
    id = "periodic_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      1 + sin(x[1])^2 + sin(x[2])^2 - 0.1 * exp(-x[1]^2 - x[2]^2)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makePeriodicFunction, "tags"),
    global.opt.params = c(0, 0),
    global.opt.value = 0.9
  )
}

class(makePeriodicFunction) = c("function", "smoof_generator")
attr(makePeriodicFunction, "name") = c("Periodic")
attr(makePeriodicFunction, "type") = c("single-objective")
attr(makePeriodicFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
