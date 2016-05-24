#' Rastrigin Function
#'
#' One of the most popular single-objective test functions consists of many
#' local optima and is thus highly multimodal with a global structure.
#' The implementation follows the formula
#' \deqn{f(\mathbf{x}) = 10n + \sum_{i=1}^{n} \left(\mathbf{x}_i^2 - 10 \cos(2\pi \mathbf{x}_i)\right).}
#' The box-constraints are given by \eqn{\mathbf{x}_i \in [-5.12, 5.12]} for
#' \eqn{i = 1, \ldots, n}.
#'
#' @references L. A. Rastrigin. Extremal control systems. Theoretical Foundations
#' of Engineering Cybernetics Series. Nauka, Moscow, 1974.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makeRastriginFunction = function(dimensions) {
  assertCount(dimensions)
  force(dimensions)
  makeSingleObjectiveFunction(
    name = paste(dimensions, "-d Rastrigin Function", sep = ""),
    id = paste0("rastrigin_", dimensions, "d"),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      n = length(x)
      10 * n + sum(x^2 - 10 * cos(2 * pi * x))
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(-5.12, dimensions),
      upper = rep(5.12, dimensions),
      vector = TRUE
    ),
    tags = attr(makeRastriginFunction, "tags"),
    global.opt.params = rep(0, dimensions),
    global.opt.value = 0
  )
}

class(makeRastriginFunction) = c("function", "smoof_generator")
attr(makeRastriginFunction, "name") = c("Rastrigin")
attr(makeRastriginFunction, "type") = c("single-objective")
attr(makeRastriginFunction, "tags") = c("single-objective", "multimodal", "continuous", "separable", "scalable")
