#' Brown Function
#'
#' This function belongs the the unimodal single-objective test functions. The
#' function is forumlated as
#' \deqn{f(\mathbf{x}) = \sum_{i = 1}^{n} (\mathbf{x}_i^2)^{(\mathbf{x}_{i + 1} + 1)} + (\mathbf{x}_{i + 1})^{(\mathbf{x}_i + 1)}}
#' subject to \eqn{\mathbf{x}_i \in [-1, 4]} for \eqn{i = 1, \ldots, n}.
#'
#' @references O. Begambre, J. E. Laier, A hybrid Particle Swarm Optimization -
#' Simplex Algorithm (PSOS) for Structural Damage Identification, Journal of
#' Advances in Engineering Software, vol. 40, no. 9, pp. 883-891, 2009.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makeBrownFunction = function(dimensions) {
  assertCount(dimensions)
  force(dimensions)
  makeSingleObjectiveFunction(
    name = "Brown Function",
    id = "brown_2d",
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      i = 1:(length(x) - 1)
      a = x[i]^2
      b = x[i + 1]^2
      sum(a^(b + 1) + b^(a + 1))
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(-1, dimensions),
      upper = rep(4, dimensions),
      vector = TRUE
    ),
    tags = attr(makeBrownFunction, "tags"),
    global.opt.params = rep(0, dimensions),
    global.opt.value = 0
  )
}

class(makeBrownFunction) = c("function", "smoof_generator")
attr(makeBrownFunction, "name") = c("Brown")
attr(makeBrownFunction, "type") = c("single-objective")
attr(makeBrownFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "scalable", "unimodal")
