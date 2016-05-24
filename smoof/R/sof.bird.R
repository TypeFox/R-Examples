#' Bird Function
#'
#' Multimodal single-objective test function. The implementation is based on the
#' mathematical formulation
#' \deqn{f(\mathbf{x}) = (\mathbf{x}_1 - \mathbf{x}_2)^2 + \exp((1 - \sin(\mathbf{x}_1))^2)\cos(\mathbf{x}_2) + \exp((1 - \cos(\mathbf{x}_2))^2)\sin(\mathbf{x}_1).}
#' The function is restricted to two dimensions with \eqn{\mathbf{x}_i \in [-2\pi, 2\pi], i = 1, 2.}
#'
#' @references S. K. Mishra, Global Optimization By Differential Evolution and
#' Particle Swarm Methods: Evaluation On Some Benchmark Functions, Munich Research
#' Papers in Economics.
#'
#' @template ret_smoof_single
#' @export
makeBirdFunction = function() {
  makeSingleObjectiveFunction(
    name = "Bird Function",
    id = "bird_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      a = (x[1] - x[2])^2
      b = exp((1 - sin(x[1]))^2) * cos(x[2])
      c = exp((1 - cos(x[2]))^2) * sin(x[1])
      return(a + b + c)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-2 * pi, -2 * pi),
      upper = c(2 * pi, 2 * pi),
      vector = TRUE
    ),
    tags = attr(makeBirdFunction, "tags"),
    global.opt.params = matrix(
      c(4.70104, 3.15294,
        -1.58214, -3.13024),
      ncol = 2, byrow = TRUE),
    global.opt.value = -106.764537
  )
}

class(makeBirdFunction) = c("function", "smoof_generator")
attr(makeBirdFunction, "name") = c("Bird")
attr(makeBirdFunction, "type") = c("single-objective")
attr(makeBirdFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
