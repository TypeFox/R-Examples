#' Price Function N. 4
#'
#' Fourth function by Price. The implementation is based on the defintion
#' \deqn{f(\mathbf{x}) = (2\mathbf{x}_1^3 \mathbf{x}_2 - \mathbf{x}_2^3)^2 + (6\mathbf{x}_1 - \mathbf{x}_2^2 + \mathbf{x}_2)^2}
#' subject to \eqn{\mathbf{x}_i \in [-500, 500]}.
#'
#' @references W. L. Price, A Controlled Random Search Procedure for Global
#' Optimisation, Computer journal, vol. 20, no. 4, pp. 367-370, 1977.
#'
#' @seealso \code{\link{makePriceN1Function}}, \code{\link{makePriceN2Function}}
#'
#' @template ret_smoof_single
#' @export
makePriceN4Function = function() {
  makeSingleObjectiveFunction(
    name = "Price Function N. 4",
    id = "price04_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      (2 * x[1]^3 * x[2] - x[2]^3)^2 + (6 * x[1] - x[2]^2 + x[2])^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-500, -500),
      upper = c(500, 500),
      vector = TRUE
    ),
    tags = attr(makePriceN4Function, "tags"),
    global.opt.params = matrix(
      c(0, 0,
        2, 4,
        1.464, -2.506),
      ncol = 2L, byrow = TRUE),
    global.opt.value = 0
  )
}

class(makePriceN4Function) = c("function", "smoof_generator")
attr(makePriceN4Function, "name") = c("Price Function N. 4")
attr(makePriceN4Function, "type") = c("single-objective")
attr(makePriceN4Function, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
