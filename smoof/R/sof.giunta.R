#' Giunta Function
#'
#' Multimodal test function based on the definition
#' \deqn{f(\mathbf{x}) = 0.6 + \sum_{i = 1}^{n} \left[\sin(\frac{16}{15} \mathbf{x}_i - 1) + \sin^2(\frac{16}{15}\mathbf{x}_i - 1) + \frac{1}{50} \sin(4(\frac{16}{15}\mathbf{x}_i - 1))\right]}
#' with box-constraints \eqn{\mathbf{x}_i \in [-1, 1]} for \eqn{i = 1, \ldots, n}.
#'
#' @references S. K. Mishra, Global Optimization By Differential Evolution and
#' Particle Swarm Methods: Evaluation On Some Benchmark Functions, Munich
#' Research Papers in Economics.
#'
#' @template ret_smoof_single
#' @export
#FIXME: this function is scalable, but global opt only known for 2D?
makeGiuntaFunction = function() {
  makeSingleObjectiveFunction(
    name = "Giunta Function",
    id = "giunta_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      a = 1.067 * x - 1
      b = sin(a)
      0.6 + sum(b + b^2 + 0.02 * sin(4 * a))
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = rep(-1, 2),
      upper = rep(1, 2),
      vector = TRUE
    ),
    tags = attr(makeGiuntaFunction, "tags"),
    global.opt.params = c(0.4673200277395354, 0.4673200169591304),
    global.opt.value = 0.06447042053690566
  )
}

class(makeGiuntaFunction) = c("function", "smoof_generator")
attr(makeGiuntaFunction, "name") = c("Giunta")
attr(makeGiuntaFunction, "type") = c("single-objective")
attr(makeGiuntaFunction, "tags") = c("single-objective", "continuous", "differentiable", "separable", "multimodal")
