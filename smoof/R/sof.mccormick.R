#' McCormick Function
#'
#' Two-dimensional, multimodal test function. The defintion is given by
#' \deqn{f(\mathbf{x}) = \sin(\mathbf{x}_1 + \mathbf{x}_2) + (\mathbf{x}_1 - \mathbf{x}_2)^2 - 1.5 \mathbf{x}_1 + 2.5 \mathbf{x}_2 + 1}
#' subject to \eqn{\mathbf{x}_1 \in [-1.5, 4], \mathbf{x}_2 \in [-3, 3]}.
#'
#' @references F. A. Lootsma (ed.), Numerical Methods for Non-Linear
#' Optimization, Academic Press, 1972.
#'
#' @template ret_smoof_single
#' @export
makeMcCormickFunction = function() {
  makeSingleObjectiveFunction(
    name = "McCormick Function",
    id = "mccormick_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      sin(x[1] + x[2]) + (x[1] - x[2])^2 - 1.5 * x[1] + 2.5 * x[2] + 1
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-1.5, -3),
      upper = c(4, 3),
      vector = TRUE
    ),
    tags = attr(makeMcCormickFunction, "tags"),
    global.opt.params = c(-0.54719, -1.54719),
    global.opt.value = -1.9133
  )
}

class(makeMcCormickFunction) = c("function", "smoof_generator")
attr(makeMcCormickFunction, "name") = c("McCormick")
attr(makeMcCormickFunction, "type") = c("single-objective")
attr(makeMcCormickFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
