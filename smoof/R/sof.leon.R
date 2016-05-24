#' Leon Function
#'
#' The function is based on the defintion
#' \deqn{f(\mathbf{x}) = 100 (\mathbf{x}_2 - \mathbf{x}_1^2)^2 + (1 - \mathbf{x}_1)^2}.
#' Box-constraints: \eqn{\mathbf{x}_i \in [-1.2, 1.2]} for \eqn{i = 1, 2}.
#'
#' @references A. Lavi, T. P. Vogel (eds), Recent Advances in Optimization
#' Techniques, John Wliley & Sons, 1966.
#'
#' @template ret_smoof_single
#' @export
makeLeonFunction = function() {
  makeSingleObjectiveFunction(
    name = "Leon Function",
    id = "leon_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-1.2, -1.2),
      upper = c(1.2, 1.2),
      vector = TRUE
    ),
    tags = attr(makeLeonFunction, "tags"),
    global.opt.params = c(1, 1),
    global.opt.value = 0
  )
}

class(makeLeonFunction) = c("function", "smoof_generator")
attr(makeLeonFunction, "name") = c("Leon")
attr(makeLeonFunction, "type") = c("single-objective")
attr(makeLeonFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "unimodal")
