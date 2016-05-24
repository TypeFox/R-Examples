#' Cube Function
#'
#' The Cube Function is defined as follows:
#' \deqn{f(\mathbf{x}) = 100 (\mathbf{x}_2 - \mathbf{x}_1^3)^2 + (1 - \mathbf{x}_1)^2.}
#' The box-constraints are given by \eqn{\mathbf{x}_i \in [-10, 10], i = 1, 2.}
#'
#' @references A. Lavi, T. P. Vogel (eds), Recent Advances in Optimization
#' Techniques, John Wliley & Sons, 1966.
#'
#' @template ret_smoof_single
#' @export
makeCubeFunction = function() {
  makeSingleObjectiveFunction(
    name = "Cube Function",
    id = "cube_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      100 * (x[2] - x[1]^3)^2 + (1 - x[1])^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makeCubeFunction, "tags"),
    global.opt.params = c(1, 1),
    global.opt.value = 0
  )
}

class(makeCubeFunction) = c("function", "smoof_generator")
attr(makeCubeFunction, "name") = c("Cube")
attr(makeCubeFunction, "type") = c("single-objective")
attr(makeCubeFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "unimodal")
