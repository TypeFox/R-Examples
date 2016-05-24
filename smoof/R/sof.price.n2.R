#' Price Function N. 2
#'
#' Second function by Price. The implementation is based on the defintion
#' \deqn{f(\mathbf{x}) = 1 + \sin^2(\mathbf{x}_1) + \sin^2(\mathbf{x}_2) - 0.1 \exp(-\mathbf{x}^2 - \mathbf{x}_2^2)}
#' subject to \eqn{\mathbf{x}_i \in [-10, 10], i = 1, 2}.
#'
#' @references W. L. Price, A Controlled Random Search Procedure for Global
#' Optimisation, Computer journal, vol. 20, no. 4, pp. 367-370, 1977.
#'
#' @seealso \code{\link{makePriceN1Function}}, \code{\link{makePriceN4Function}}
#'
#' @template ret_smoof_single
#' @export
makePriceN2Function = function() {
  makeSingleObjectiveFunction(
    name = "Price Function N. 2",
    id = "price02_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      1 + sin(x[1])^2 + sin(x[2])^2 - 0.1 * exp(-x[1]^2 - x[2]^2)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makePriceN2Function, "tags"),
    global.opt.params = c(0, 0),
    global.opt.value = 0.9
  )
}

class(makePriceN2Function) = c("function", "smoof_generator")
attr(makePriceN2Function, "name") = c("Price N. 2")
attr(makePriceN2Function, "type") = c("single-objective")
attr(makePriceN2Function, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
