#' Bukin function N. 6
#'
#' Beside Bukin N. 2 and N. 4 this is the last \dQuote{Bukin family} function.
#' It is given by the formula
#' \deqn{f(\mathbf{x}) = 100 \sqrt{||\mathbf{x}_2 - 0.01 \mathbf{x}_1^2||} + 0.01 ||\mathbf{x}_1 + 10||}
#' and the box constraints \eqn{mathbf{x}_1 \in [-15, 5], \mathbf{x}_2 \in [-3, 3]}.
#'
#' @references Z. K. Silagadze, Finding Two-Dimesnional Peaks, Physics of Particles
#' and Nuclei Letters, vol. 4, no. 1, pp. 73-80, 2007.
#'
#' @seealso \code{\link{makeBukinN2Function}}, \code{\link{makeBukinN4Function}}
#'
#' @template ret_smoof_single
#' @export
makeBukinN6Function = function() {
  makeSingleObjectiveFunction(
    name = "Bukin Function N.6",
    id = "bukin06_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      100 * sqrt(abs(x[2] - 0.01 * x[1]^2)) + 0.01 * abs(x[1] + 10)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-15, -3),
      upper = c(-5, 3),
      vector = TRUE
    ),
    tags = attr(makeBukinN6Function, "tags"),
    global.opt.params = c(-10, 1),
    global.opt.value = 0
  )
}

class(makeBukinN6Function) = c("function", "smoof_generator")
attr(makeBukinN6Function, "name") = c("Bukin N. 6")
attr(makeBukinN6Function, "type") = c("single-objective")
attr(makeBukinN6Function, "tags") = c("single-objective", "continuous", "non-differentiable", "non-separable", "non-scalable", "multimodal")
