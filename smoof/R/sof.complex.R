#' @title Complex function.
#'
#' @description Two-dimensional test function based on the formula
#' \deqn{f(\mathbf{x}) = (x_1^3 - 3 x_1 x_2^2 - 1)^2 + (3 x_2 x_1^2 - x_2^3)^2}
#' with \eqn{\mathbf{x}_1, \mathbf{x}_2 \in [-2, 2]}.
#'
#' @references See \url{http://al-roomi.org/benchmarks/unconstrained/2-dimensions/43-complex-function}.
#'
#' @template ret_smoof_single
#' @export
makeComplexFunction = function() {
  makeSingleObjectiveFunction(
    name = "Complex Function",
    id = paste0("complex_2d"),
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      (x[1]^3 - 3 * x[1] * x[2]^2 - 1)^2 + (3 * x[2] * x[1]^2 - x[2]^3)^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-2, -2),
      upper = c(2, 2),
      vector = TRUE
    ),
    tags = attr(makeComplexFunction, "tags"),
    global.opt.params = c(1, 0),
    global.opt.value = 0
  )
}

class(makeComplexFunction) = c("function", "smoof_generator")
attr(makeComplexFunction, "name") = c("Complex")
attr(makeComplexFunction, "type") = c("single-objective")
attr(makeComplexFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
