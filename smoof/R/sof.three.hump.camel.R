#' Three-Hump Camel Function
#'
#' This two-dimensional function is based on the defintion
#' \deqn{f(\mathbf{x}) = 2 \mathbf{x}_1^2 - 1.05 \mathbf{x}_1^4 + \frac{\mathbf{x}_1^6}{6} + \mathbf{x}_1\mathbf{x}_2 + \mathbf{x}_2^2}
#' subject to \eqn{-5 \leq \mathbf{x}_i \leq 5}.
#'
#' @references F. H. Branin Jr., Widely Convergent Method of Finding Multiple
#' Solutions of Simul- taneous Nonlinear Equations, IBM Journal of Research
#' and Development, vol. 16, no. 5, pp. 504-522, 1972.
#'
#' @template ret_smoof_single
#' @export
makeThreeHumpCamelFunction = function() {
  makeSingleObjectiveFunction(
    name = "Three-Hump Camel Function",
    id = "threeHumpCamel_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      2 * x[1]^2 - 1.05 * x[1]^4 + (x[1]^6) / 6 + x[1] * x[2] + x[2]^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-5, -5),
      upper = c(5, 5),
      vector = TRUE
    ),
    tags = attr(makeThreeHumpCamelFunction, "tags"),
    global.opt.params = c(0, 0),
    global.opt.value = 0
  )
}

class(makeThreeHumpCamelFunction) = c("function", "smoof_generator")
attr(makeThreeHumpCamelFunction, "name") = c("Three-Hump Camel")
attr(makeThreeHumpCamelFunction, "type") = c("single-objective")
attr(makeThreeHumpCamelFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
