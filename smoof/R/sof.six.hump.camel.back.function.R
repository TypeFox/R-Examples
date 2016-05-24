#' Three-Hump Camel Function
#'
#' Two dimensional single-objective test function with six local minima oh which
#' two are global. The surface is similar to the back of a camel. That is why it
#' is called Camel function. The implementation is based on the formula:
#' \deqn{f(\mathbf{x}) = \left(4 - 2.1\mathbf{x}_1^2 + \mathbf{x}_1^{0.75}\right)\mathbf{x}_1^2 + \mathbf{x}_1 \mathbf{x}_2 + \left(-4 + 4\mathbf{x}_2^2\right)\mathbf{x}_2^2}
#' with box constraints \eqn{\mathbf{x}_1 \in [-3, 3]} and \eqn{\mathbf{x}_2 \in [-2, 2]}.
#'
#' @references Dixon, L. C. W. and Szego, G. P.: The optimization problem: An introduction.
#' In: Towards Global Optimization II, New York: North Holland, 1978.
#'
#' @template ret_smoof_single
#' @export
makeSixHumpCamelFunction = function() {
  makeSingleObjectiveFunction(
    name = "Six-Hump Camel Back Function",
    id = "sixHumpCamelBack_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      xx1 = x[1]^2
      xx2 = x[2]^2
      (4 - 2.1 * xx1 + xx1^2 / 3) * xx1 + x[1] * x[2] + (4 * xx2 - 4) * xx2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-3, -2),
      upper = c(3, 2),
      vector = TRUE
    ),
    tags = attr(makeSixHumpCamelFunction, "tags"),
    global.opt.params = matrix(
      c(-0.0898, 0.7126,
        0.0898, -0.7126),
      ncol = 2L, byrow = TRUE),
    global.opt.value = -1.0316
  )
}

class(makeSixHumpCamelFunction) = c("function", "smoof_generator")
attr(makeSixHumpCamelFunction, "name") = c("Six-Hump Camel Back")
attr(makeSixHumpCamelFunction, "type") = c("single-objective")
attr(makeSixHumpCamelFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
