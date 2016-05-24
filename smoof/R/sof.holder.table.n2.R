#' Holder Table function N. 2
#'
#' This multimodal function is defined as
#' \deqn{f(\mathbf{x}) = -\left|\sin(\mathbf{x}_1)\cos(\mathbf{x}_2)\exp(|1 - \sqrt{\mathbf{x}_1 + \mathbf{x}_2}/\pi|)\right|}
#' with box-constraints \eqn{\mathbf{x}_i \in [-10, 10], i = 1, 2}.
#'
#' @references S. K. Mishra, Global Optimization By Differential Evolution and
#' Particle Swarm Methods: Evaluation On Some Benchmark Functions, Munich
#' Research Papers in Economics.
#'
#' @seealso \code{\link{makeHolderTableN1Function}}
#'
#' @template ret_smoof_single
#' @export
makeHolderTableN2Function = function() {
  makeSingleObjectiveFunction(
    name = "Holder Table Function N. 2",
    id = "holderTable02_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      -abs(sin(x[1]) * cos(x[2]) * exp(abs(1 - sqrt(x[1]^2 + x[2]^2) / 3.1415)))
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makeHolderTableN2Function, "tags"),
    global.opt.params = matrix(
      c(8.05502, 9.66459,
        -8.05502, 9.66459,
        8.05502, -9.66459,
        -8.05502, -9.66459),
      ncol = 2L, byrow = TRUE),
    global.opt.value = -19.2085
  )
}

class(makeHolderTableN2Function) = c("function", "smoof_generator")
attr(makeHolderTableN2Function, "name") = c("Holder Table N. 2")
attr(makeHolderTableN2Function, "type") = c("single-objective")
attr(makeHolderTableN2Function, "tags") = c("single-objective", "continuous", "differentiable", "separable", "non-scalable", "multimodal")
