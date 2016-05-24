#' Adjiman function
#'
#' This two-dimensional multimodal test function follows the formula
#' \deqn{f(\mathbf{x}) = \cos(\mathbf{x}_1)\sin(\mathbf{x}_2) - \frac{\mathbf{x}_1}{(\mathbf{x}_2^2 + 1)}}
#' with \eqn{\mathbf{x}_1 \in [-1, 2], \mathbf{x}_2 \in [2, 1]}.
#'
#' @references C. S. Adjiman, S. Sallwig, C. A. Flouda, A. Neumaier, A Global Optimization
#' Method, aBB for General Twice-Differentiable NLPs-1, Theoretical Advances, Computers
#' Chemical Engineering, vol. 22, no. 9, pp. 1137-1158, 1998.
#'
#' @template ret_smoof_single
#' @export
makeAdjimanFunction = function() {
  makeSingleObjectiveFunction(
    name = "Adjiman Function",
    id = "adjiman_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      cos(x[1]) * sin(x[2]) - x[1] / (x[2]^2 + 1)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-1, -1),
      upper = c(2, 1),
      vector = TRUE
    ),
    tags = attr(makeAdjimanFunction, "tags"),
    global.opt.params = c(2, 0.10578),
    global.opt.value = -2.02181
  )
}

class(makeAdjimanFunction) = c("function", "smoof_generator")
attr(makeAdjimanFunction, "name") = c("Adjiman")
attr(makeAdjimanFunction, "type") = c("single-objective")
attr(makeAdjimanFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
