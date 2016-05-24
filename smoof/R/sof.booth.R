#' Booth Function
#'
#' This function is based on the formula
#' \deqn{f(\mathbf{x}) = (\mathbf{x}_1 + 2\mathbf{x}_2 - 7)^2 + (2\mathbf{x}_1 + \mathbf{x}_2 - 5)^2}
#' subject to \eqn{\mathbf{x}_i \in [-10, 10], i = 1, 2}.
#'
#' @template ret_smoof_single
#' @export
makeBoothFunction = function() {
  makeSingleObjectiveFunction(
    name = "Booth Function",
    id = "booth_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      (x[1] + 2 * x[2] - 7)^2 + (2 * x[1] + x[2] - 5)^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makeBoothFunction, "tags"),
    global.opt.params = c(1, 3),
    global.opt.value = 0
  )
}

class(makeBoothFunction) = c("function", "smoof_generator")
attr(makeBoothFunction, "name") = c("Booth")
attr(makeBoothFunction, "type") = c("single-objective")
attr(makeBoothFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "unimodal")
