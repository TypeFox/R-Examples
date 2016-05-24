#' Generalized Drop-Wave Function
#'
#' Multimodal single-objective function following the formula:
#' \deqn{\mathbf{x} = -\frac{1 + \cos(\sqrt{\sum_{i = 1}^{n} \mathbf{x}_i^2})}{2 + 0.5 \sum_{i = 1}^{n} \mathbf{x}_i^2}}
#' with \eqn{\mathbf{x}_i \in [-5.12, 5.12], i = 1, \ldots, n}.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makeGeneralizedDropWaveFunction = function(dimensions = 2L) {
  assertCount(dimensions)
  force(dimensions)
  makeSingleObjectiveFunction(
    name = paste(dimensions, "-d Generelized Drop-Wave Function", sep = ""),
    id = paste0("generalizedDropWave_", dimensions, "d"),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      a = sum(x^2)
      -(1 + cos(12 * sqrt(a))) / (0.5 * a + 2)
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(-5.12, dimensions),
      upper = rep(5.12, dimensions),
      vector = TRUE
    ),
    tags = attr(makeGeneralizedDropWaveFunction, "tags"),
    global.opt.params = rep(0, dimensions),
    global.opt.value = -1
  )
}

class(makeGeneralizedDropWaveFunction) = c("function", "smoof_generator")
attr(makeGeneralizedDropWaveFunction, "name") = c("Generelized Drop-Wave")
attr(makeGeneralizedDropWaveFunction, "type") = c("single-objective")
attr(makeGeneralizedDropWaveFunction, "tags") = c("single-objective", "multimodal", "non-separable", "continuous", "differentiable", "scalable")
