#' Bohachevsky function N. 1
#'
#' Highly multimodal single-objective test function. The mathematical formula is
#' given by
#' \deqn{f(\mathbf{x}) = \sum_{i = 1}^{n - 1} (\mathbf{x}_i^2 + 2 \mathbf{x}_{i + 1}^2 - 0.3\cos(3\pi\mathbf{x}_i) - 0.4\cos(4\pi\mathbf{x}_{i + 1}) + 0.7)}
#' with box-constraints \eqn{\mathbf{x}_i  \in [-100, 100]} for \eqn{i = 1, \ldots, n}.
#' The multimodality will be visible by \dQuote{zooming in} in the plot.
#'
#' @references  I. O. Bohachevsky, M. E. Johnson, M. L. Stein, General Simulated
#' Annealing for Function Optimization, Technometrics, vol. 28, no. 3, pp. 209-217, 1986.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makeBohachevskyN1Function = function(dimensions) {
  assertCount(dimensions)
  force(dimensions)
  makeSingleObjectiveFunction(
    name = paste(dimensions, "-d Bohachevsky Function N. 1", sep = ""),
    id = paste0("bohachevsky01_", dimensions, "d"),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      i = 1:(length(x) - 1)
      sum(x[i]^2 + 2 * x[i + 1]^2 - 0.3 * cos(3 * pi * x[i]) - 0.4 * cos(4 * pi * x[i + 1]) + 0.7)
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(-15, dimensions),
      upper = rep(15, dimensions),
      vector = TRUE
    ),
    tags = attr(makeBohachevskyN1Function, "tags"),
    global.opt.params = rep(0, dimensions),
    global.opt.value = 0
  )
}

class(makeBohachevskyN1Function) = c("function", "smoof_generator")
attr(makeBohachevskyN1Function, "name") = c("Bohachevsky N. 1")
attr(makeBohachevskyN1Function, "type") = c("single-objective")
attr(makeBohachevskyN1Function, "tags") = c("single-objective", "continuous", "differentiable", "separable", "scalable", "multimodal")
