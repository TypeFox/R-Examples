#' Cosine Mixture Function
#'
#' Single-objective test function based on the formula
#' \deqn{f(\mathbf{x}) = -0.1 \sum_{i = 1}^{n} \cos(5\pi\mathbf{x}_i) - \sum_{i = 1}^{n} \mathbf{x}_i^2}
#' subject to \eqn{\mathbf{x}_i \in [-1, 1]} for \eqn{i = 1, \ldots, n}.
#'
#' @references M. M. Ali, C. Khompatraporn, Z. B. Zabinsky, A Numerical Evaluation
#' of Several Stochastic Algorithms on Selected Continuous Global Optimization
#' Test Problems, Journal of Global Optimization, vol. 31, pp. 635-672, 2005.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makeCosineMixtureFunction = function(dimensions) {
  assertCount(dimensions)
  force(dimensions)
  makeSingleObjectiveFunction(
    name = "Cosine Mixture Function",
    id = paste0("cosineMixture", dimensions, "d"),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      a = -0.1 * sum(cos(5 * pi * x))
      b = sum(x^2)
      return (a - b)
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(-1, dimensions),
      upper = rep(1, dimensions),
      vector = TRUE
    ),
    minimize = FALSE,
    tags = attr(makeCosineMixtureFunction, "tags"),
    global.opt.params = rep(0, dimensions),
    global.opt.value = -0.1 * dimensions
  )
}

class(makeCosineMixtureFunction) = c("function", "smoof_generator")
attr(makeCosineMixtureFunction, "name") = c("Cosine Mixture")
attr(makeCosineMixtureFunction, "type") = c("single-objective")
attr(makeCosineMixtureFunction, "tags") = c("single-objective", "discontinuous", "non-differentiable", "separable", "scalable", "multimodal")
