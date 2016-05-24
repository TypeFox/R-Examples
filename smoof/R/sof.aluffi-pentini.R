#' @title Aluffi-Pentini function.
#'
#' @description Two-dimensional test function based on the formula
#' \deqn{f(\mathbf{x}) = 0.25 x_1^4 - 0.5 x_1^2 + 0.1 x_1 + 0.5 x_2^2}
#' with \eqn{\mathbf{x}_1, \mathbf{x}_2 \in [-10, 10]}.
#'
#' @references See \url{http://al-roomi.org/benchmarks/unconstrained/2-dimensions/26-aluffi-pentini-s-or-zirilli-s-function}.
#'
#' @note This functions is also know as the Zirilli function.
#'
#' @template ret_smoof_single
#' @export
makeAluffiPentiniFunction = function() {
  makeSingleObjectiveFunction(
    name = "Aluffi-Pentini Function",
    id = "allufiPentini_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      0.25 * x[1]^4 - 0.5 * x[1]^2 + 0.1 * x[1] + 0.5 * x[2]^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makeAluffiPentiniFunction, "tags"),
    global.opt.params = c(-1.046680576580755, 0),
    global.opt.value = -0.352386073800034
  )
}

class(makeAluffiPentiniFunction) = c("function", "smoof_generator")
attr(makeAluffiPentiniFunction, "name") = c("Aluffi-Pentini")
attr(makeAluffiPentiniFunction, "type") = c("single-objective")
attr(makeAluffiPentiniFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "unimodal")
