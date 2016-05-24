#' @title Kearfott function.
#'
#' @description Two-dimensional test function based on the formula
#' \deqn{f(\mathbf{x}) = (x_1^2 + x_2^2 - 2)^2 + (x_1^2 - x_2^2 - 1)^2}
#' with \eqn{\mathbf{x}_1, \mathbf{x}_2 \in [-3, 4]}.
#'
#' @references See \url{http://al-roomi.org/benchmarks/unconstrained/2-dimensions/59-kearfott-s-function}.
#'
#' @template ret_smoof_single
#' @export
makeKearfottFunction = function() {
  makeSingleObjectiveFunction(
    name = "Kearfott Function",
    id = "kearfott_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      (x[1]^2 + x[2]^2 - 2)^2 + (x[1]^2 - x[2]^2 - 1)^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-3, -3),
      upper = c(4, 4),
      vector = TRUE
    ),
    tags = attr(makeKearfottFunction, "tags"),
    global.opt.params = matrix(
      c(1.22474487139159, 0.707106781186548,
        -1.22474487139159, -0.707106781186548,
        -1.22474487139159, 0.707106781186548,
        1.22474487139159, -0.707106781186548),
      ncol = 2, byrow = TRUE),
    global.opt.value = 0
  )
}

class(makeKearfottFunction) = c("function", "smoof_generator")
attr(makeKearfottFunction, "name") = c("Kearfott")
attr(makeKearfottFunction, "type") = c("single-objective")
attr(makeKearfottFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
