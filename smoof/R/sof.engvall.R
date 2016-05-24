#' @title Complex function.
#'
#' @description Two-dimensional test function based on the formula
#' \deqn{f(\mathbf{x}) = (x_1^4 + x_2^4 + 2 x_1^2 x_2^2 - 4 x_1 + 3}
#' with \eqn{\mathbf{x}_1, \mathbf{x}_2 \in [-2000, 2000]}.
#'
#' @references See \url{http://al-roomi.org/benchmarks/unconstrained/2-dimensions/116-engvall-s-function}.
#'
#' @template ret_smoof_single
#' @export
makeEngvallFunction = function() {
  makeSingleObjectiveFunction(
    name = "Engvall Function",
    id = "engvall_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      x[1]^4 + x[2]^4 + 2 * x[1]^2 * x[2]^2 - 4 * x[1] + 3
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-2000, -2000),
      upper = c(2000, 2000),
      vector = TRUE
    ),
    tags = attr(makeEngvallFunction, "tags"),
    global.opt.params = c(1, 0),
    global.opt.value = 0
  )
}

class(makeEngvallFunction) = c("function", "smoof_generator")
attr(makeEngvallFunction, "name") = c("Engvall")
attr(makeEngvallFunction, "type") = c("single-objective")
attr(makeEngvallFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "unimodal")
