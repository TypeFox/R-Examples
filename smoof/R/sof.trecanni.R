#' Trecanni Function
#'
#' The Trecanni function belongs to the unimodal test functions. It is based
#' on the formula
#' \deqn{f(\mathbf(x)) = \mathbf(x)_1^4 - 4 \mathbf(x)_1^3 + 4 \mathbf(x)_1 + \mathbf(x)_2^2.}
#' The box-constraints \eqn{\mathbf{x}_i \in [-5, 5], i = 1, 2} define the
#' domain of defintion.
#'
#' @references L. C. W. Dixon, G. P. Szego (eds.), Towards Global Optimization 2,
#' Elsevier, 1978.
#'
#' @template ret_smoof_single
#' @export
makeTrecanniFunction = function() {
  makeSingleObjectiveFunction(
    name = "Trecanni Function",
    id = "trecanni_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      x[1]^4 + 4 * (x[1]^3 + x[1]^2) + x[2]^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-5, -5),
      upper = c(5, 5),
      vector = TRUE
    ),
    tags = c("continuous", "differentiable", "separable", "non-scalable", "unimodal"),
    global.opt.params = matrix(
      c(0, 0,
        -2, 0),
      ncol = 2L, byrow = TRUE),
    global.opt.value = 0
  )
}

class(makeTrecanniFunction) = c("function", "smoof_generator")
attr(makeTrecanniFunction, "name") = c("Trecanni")
attr(makeTrecanniFunction, "type") = c("single-objective")
attr(makeTrecanniFunction, "tags") = c("single-objective", "continuous", "differentiable", "separable", "non-scalable", "unimodal")
