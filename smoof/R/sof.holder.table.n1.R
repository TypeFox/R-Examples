#' Holder Table function N. 1
#'
#' This multimodal function is defined as
#' \deqn{f(\mathbf{x}) = -\left|\cos(\mathbf{x}_1)\cos(\mathbf{x}_2)\exp(|1 - \sqrt{\mathbf{x}_1 + \mathbf{x}_2}/\pi|)\right|}
#' with box-constraints \eqn{\mathbf{x}_i \in [-10, 10], i = 1, 2}.
#'
#' @references S. K. Mishra, Global Optimization By Differential Evolution and
#' Particle Swarm Methods: Evaluation On Some Benchmark Functions, Munich
#' Research Papers in Economics.
#'
#' @seealso \code{\link{makeHolderTableN2Function}}
#'
#' @template ret_smoof_single
#' @export
makeHolderTableN1Function = function() {
  makeSingleObjectiveFunction(
    name = "Holder Table Function N. 1",
    id = "holderTable01_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      -abs(cos(x[1]) * cos(x[2]) * exp(abs(1 - sqrt(x[1]^2 + x[2]^2) / 3.1415)))
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makeHolderTableN1Function, "tags"),
    global.opt.params = matrix(
      c(9.646168, 9.646168,
        -9.646168, 9.646168,
        9.646168, -9.646168,
        -9.646168, -9.646168),
      ncol = 2L, byrow = TRUE),
    global.opt.value = -26.920336
  )
}

class(makeHolderTableN1Function) = c("function", "smoof_generator")
attr(makeHolderTableN1Function, "name") = c("Holder Table N. 1")
attr(makeHolderTableN1Function, "type") = c("single-objective")
attr(makeHolderTableN1Function, "tags") = c("single-objective", "continuous", "differentiable", "separable", "non-scalable", "multimodal")
