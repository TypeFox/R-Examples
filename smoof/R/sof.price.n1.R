#' Price Function N. 1
#'
#' Second function by Price. The implementation is based on the defintion
#' \deqn{f(\mathbf{x}) = (|\mathbf{x}_1| - 5)^2 + (|\mathbf{x}_2 - 5)^2}
#' subject to \eqn{\mathbf{x}_i \in [-500, 500], i = 1, 2}.
#'
#' @references W. L. Price, A Controlled Random Search Procedure for Global
#' Optimisation, Computer journal, vol. 20, no. 4, pp. 367-370, 1977.
#'
#' @seealso \code{\link{makePriceN2Function}}, \code{\link{makePriceN4Function}}
#'
#' @template ret_smoof_single
#' @export
makePriceN1Function = function() {
  makeSingleObjectiveFunction(
    name = "Price Function N. 1",
    id = "price01_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      sum((abs(x) - 5)^2)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-500, -500),
      upper = c(500, 500),
      vector = TRUE
    ),
    tags = attr(makePriceN1Function, "tags"),
    global.opt.params = matrix(
      c(5, 5,
        -5, 5,
        5, -5,
        -5, 5),
      ncol = 2L, byrow = TRUE),
    global.opt.value = 0
  )
}

class(makePriceN1Function) = c("function", "smoof_generator")
attr(makePriceN1Function, "name") = c("Price N. 1")
attr(makePriceN1Function, "type") = c("single-objective")
attr(makePriceN1Function, "tags") = c("single-objective", "continuous", "non-differentiable", "separable", "non-scalable", "multimodal")
