#' Branin RCOS function
#'
#' Popular 2-dimensional single-objective test function based on the formula:
#' \deqn{f(\mathbf{x}) = a \left(\mathbf{x}_2 - b \mathbf{x}_1^2 + c \mathbf{x_1} - d\right)^2 + e\left(1 - f\right)\cos(\mathbf{x}_1) + e,}
#' where \eqn{a = 1, b = \frac{5.1}{4\pi^2}, c = \frac{5}{\pi}, d = 6, e = 10} and
#' \eqn{f = \frac{1}{8\pi}}. The box constraints are given by \eqn{\mathbf{x}_1 \in [-5, 10]}
#' and \eqn{\mathbf{x}_2 \in [0, 15]}. The function has three global minima.
#'
#' @references F. H. Branin. Widely convergent method for finding
#' multiple solutions of simultaneous nonlinear equations.
#' IBM J. Res. Dev. 16, 504-522, 1972.
#'
#' @examples
#' library(ggplot2)
#' fn = makeBraninFunction()
#' print(fn)
#' print(autoplot(fn, show.optimum = TRUE))
#' @template ret_smoof_single
#' @export
makeBraninFunction = function() {
  makeSingleObjectiveFunction(
    name = "Branin RCOS Function",
    id = "branin_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      a = 1
      b = 5.1 / (4 * pi^2)
      c = 5 / pi
      d = 6
      e = 10
      f = 1 / (8 * pi)
      return (a * (x[2] - b * x[1]^2 + c * x[1] - d)^2 + e * (1 - f) * cos(x[1]) + e)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-5, 0),
      upper = c(10, 15),
      vector = TRUE
    ),
    tags = attr(makeBraninFunction, "tags"),
    global.opt.params = matrix(
      c(-pi, 12.275, pi, 2.275, 3 * pi, 2.475),
      ncol = 2L, byrow = TRUE),
    global.opt.value = 0.397887
  )
}

class(makeBraninFunction) = c("function", "smoof_generator")
attr(makeBraninFunction, "name") = c("BraninRCOS")
attr(makeBraninFunction, "type") = c("single-objective")
attr(makeBraninFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
