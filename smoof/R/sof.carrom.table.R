#' Carrom Table Function
#'
#' This function is defined as follows:
#' \deqn{f(\mathbf{x}) = -\frac{1}{30} \left((\cos(\mathbf{x}_1)\exp(|1 - ((\mathbf{x}_1^2 + \mathbf{x}_2^2)^{0.5} / \pi)^2)|\right).}
#' The box-constraints are given by \eqn{\mathbf{x}_i \in [-10, 10], i = 1, 2}.
#'
#' @references S. K. Mishra, Global Optimization By Differential Evolution and
#' Particle Swarm Methods: Evaluation On Some Benchmark Functions, Munich
#' Research Papers in Economics.
#'
#' @template ret_smoof_single
#' @export
makeCarromTableFunction = function() {
  makeSingleObjectiveFunction(
    name = "Carrom Table Function",
    id = "carromTable_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      (-1 / 30) * exp(2 * abs(1 - (sqrt(x[1]^2 + x[2]^2) / pi))) * cos(x[1])^2 * cos(x[2])^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makeCarromTableFunction, "tags"),
    global.opt.params = matrix(
      c(9.646157266348881, 9.646134286497169,
        -9.646157266348881, 9.646134286497169,
        9.646157266348881, -9.646134286497169,
        -9.646157266348881, -9.646134286497169),
      ncol = 2L, byrow = TRUE),
    global.opt.value = -24.1568155
  )
}

class(makeCarromTableFunction) = c("function", "smoof_generator")
attr(makeCarromTableFunction, "name") = c("Carrom Table")
attr(makeCarromTableFunction, "type") = c("single-objective")
attr(makeCarromTableFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
