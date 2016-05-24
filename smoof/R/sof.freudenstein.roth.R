#' Freudenstein Roth Function
#'
#' This test function is based on the formula
#' \deqn{f(\mathbf{x}) = (\mathbf{x}_1 - 13 + ((5 - \mathbf{x}_2)\mathbf{x}_2 - 2)\mathbf{x}_2)^2 + (\mathbf{x}_1 - 29 + ((\mathbf{x}_2 + 1)\mathbf{x}_2 - 14)\mathbf{x}_2)^2}
#' subject to \eqn{\mathbf{x}_i \in [-10, 10], i = 1, 2}.
#'
#' @references S. S. Rao, Engineering Optimization: Theory and Practice,
#' John Wiley & Sons, 2009.
#'
#' @template ret_smoof_single
#' @export
makeFreudensteinRothFunction = function() {
  makeSingleObjectiveFunction(
    name = "Freudenstein Roth Function",
    id = "freudensteinRoth_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      (x[1] - 13 + ((5 - x[2]) * x[2] - 2) * x[2])^2 + (x[1] - 29 + ((x[2] + 1) * x[2] - 14) * x[2])^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makeFreudensteinRothFunction, "tags"),
    global.opt.params = c(5, 4),
    global.opt.value = 0
  )
}

class(makeFreudensteinRothFunction) = c("function", "smoof_generator")
attr(makeFreudensteinRothFunction, "name") = c("Freudenstein Roth")
attr(makeFreudensteinRothFunction, "type") = c("single-objective")
attr(makeFreudensteinRothFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
