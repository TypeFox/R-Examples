#' Zettl Function
#'
#' The unimodal Zettl Function is based on the defintion
#' \deqn{f(\mathbf{x}) = (\mathbf{x}_1^2 + \mathbf{x}_2^2 - 2\mathbf{x}_1)^2 + 0.25 \mathbf{x}_1}
#' with box-constraints \eqn{\mathbf{x}_i \in [-5, 10], i = 1, 2}.
#'
#' @references H. P. Schwefel, Evolution and Optimum Seeking, John Wiley Sons, 1995.
#'
#' @template ret_smoof_single
#' @export
makeZettlFunction = function() {
  makeSingleObjectiveFunction(
    name = "Zettl Function",
    id = "zettl_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      (x[1]^2 + x[2]^2 - 2 * x[1])^2 + 0.25 * x[1]
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-5, -5),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = c("continuous", "differentiable", "non-separable", "non-scalable", "unimodal"),
    global.opt.params = c(-0.0299, 0),
    global.opt.value = -0.003791
  )
}

class(makeZettlFunction) = c("function", "smoof_generator")
attr(makeZettlFunction, "name") = c("Zettl")
attr(makeZettlFunction, "type") = c("single-objective")
attr(makeZettlFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "unimodal")
