#' Bukin function N. 4
#'
#' Second continous Bukin function test function. The formula is given by
#' \deqn{f(\mathbf{x}) = 100 \mathbf{x}_2^2 + 0.01 * ||\mathbf{x}_1 +10||}
#' and the box constraints \eqn{mathbf{x}_1 \in [-15, 5], \mathbf{x}_2 \in [-3, 3]}.
#'
#' @references Z. K. Silagadze, Finding Two-Dimesnional Peaks, Physics of Particles
#' and Nuclei Letters, vol. 4, no. 1, pp. 73-80, 2007.
#'
#' @seealso \code{\link{makeBukinN2Function}}, \code{\link{makeBukinN6Function}}
#'
#' @template ret_smoof_single
#' @export
makeBukinN4Function = function() {
  makeSingleObjectiveFunction(
    name = "Bukin Function N. 4",
    id = "bukin04_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      100 * x[2]^2 + 0.01 * abs(x[1] + 10)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-15, -3),
      upper = c(-5, 3),
      vector = TRUE
    ),
    tags = attr(makeBukinN4Function, "tags"),
    global.opt.params = c(-10, 0),
    global.opt.value = 0
  )
}

class(makeBukinN4Function) = c("function", "smoof_generator")
attr(makeBukinN4Function, "name") = c("Bukin N. 4")
attr(makeBukinN4Function, "type") = c("single-objective")
attr(makeBukinN4Function, "tags") = c("single-objective", "continuous", "non-differentiable", "separable", "non-scalable", "multimodal")
