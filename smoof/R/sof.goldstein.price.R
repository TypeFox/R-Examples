#' Goldstein-Price Function
#'
#' Two-dimensional test function for global optimization. The implementation
#' follows the formula:
#' \deqn{f(\mathbf{x}) = \left(1 + (\mathbf{x}_1 + \mathbf{x}_2 + 1)^2 \cdot (19 - 14\mathbf{x}_1 + 3\mathbf{x}_1^2 - 14\mathbf{x}_2 + 6\mathbf{x}_1\mathbf{x}_2 + 3\mathbf{x}_2^2)\right)\\ \qquad \cdot \left(30 + (2\mathbf{x}_1 - 3\mathbf{x}_2)^2 \cdot (18 - 32\mathbf{x}_1 + 12\mathbf{x}_1^2 + 48\mathbf{x}_2 - 36\mathbf{x}_1\mathbf{x}_2 + 27\mathbf{x}_2^2)\right)}
#' with \eqn{\mathbf{x}_i \in [-2, 2], i = 1, 2}.
#'
#' @references Goldstein, A. A. and Price, I. F.: On descent from local minima.
#' Math. Comput., Vol. 25, No. 115, 1971.
#'
#' @template ret_smoof_single
#' @export
makeGoldsteinPriceFunction = function() {
  makeSingleObjectiveFunction(
    name = "Goldstein-Price Function",
    id = "goldsteinPrice_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      xx1 = x[1]^2
      xx2 = x[2]^2
      xx12 = x[1] * x[2]
      a = 1 + (x[1] + x[2] + 1)^2 * (19 - 14 * x[1] + 3 * xx1 - 14 * x[2] + 6 * xx12 + 3 * xx2)
      b = 30 + (2 * x[1] - 3 * x[2])^2 * (18 - 32 * x[1] + 12 * xx1 + 48 * x[2] - 36 * xx12 + 27 * xx2)
      return (a * b)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-2, -2),
      upper = c(2, 2),
      vector = TRUE
    ),
    tags = attr(makeGoldsteinPriceFunction, "tags"),
    global.opt.params = c(0, -1),
    global.opt.value = 3
  )
}

class(makeGoldsteinPriceFunction) = c("function", "smoof_generator")
attr(makeGoldsteinPriceFunction, "name") = c("Goldstein-Price")
attr(makeGoldsteinPriceFunction, "type") = c("single-objective")
attr(makeGoldsteinPriceFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
