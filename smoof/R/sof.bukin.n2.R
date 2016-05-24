#' Bukin function N. 2
#'
#' Muldimodal, non-scalable, continous optimization test function given by:
#' \deqn{f(\mathbf{x}) = 100 (\mathbf{x}_2 - 0.01 * \mathbf{x}_1^2 + 1) + 0.01 (\mathbf{x}_1 + 10)^2}
#' subject to \eqn{\mathbf{x}_1 \in [-15, 5]} and \eqn{\mathbf{x}_2 \in [-3, 3]}.
#'
#' @references Z. K. Silagadze, Finding Two-Dimesnional Peaks, Physics of Particles
#' and Nuclei Letters, vol. 4, no. 1, pp. 73-80, 2007.
#'
#' @seealso \code{\link{makeBukinN4Function}}, \code{\link{makeBukinN6Function}}
#'
#' @template ret_smoof_single
#' @export
makeBukinN2Function = function() {
  makeSingleObjectiveFunction(
    name = "Bukin Function N. 2",
    id = "bukin02_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      100 * (x[2] - 0.01 * x[1]^2 + 1) + 0.01 * (x[1] + 10)^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-15, -3),
      upper = c(-5, 3),
      vector = TRUE
    ),
    tags = attr(makeBukinN2Function, "tags"),
    global.opt.params = c(-10, 0),
    global.opt.value = 0
  )
}

class(makeBukinN2Function) = c("function", "smoof_generator")
attr(makeBukinN2Function, "name") = c("Bukin N. 2")
attr(makeBukinN2Function, "type") = c("single-objective")
attr(makeBukinN2Function, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
