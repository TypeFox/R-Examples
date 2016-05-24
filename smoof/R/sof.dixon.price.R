#' Dixon-Price Function
#'
#' Dixon and Price defined the function
#' \deqn{f(\mathbf{x}) = (\mathbf{x}_1 - 1)^2 + \sum_{i = 1}^{n} i (2\mathbf{x}_i^2 - \mathbf{x}_{i - 1})}
#' subject to \eqn{\mathbf{x}_i \in [-10, 10]} for \eqn{i = 1, \ldots, n}.
#'
#' @references L. C. W. Dixon, R. C. Price, The Truncated Newton Method for
#' Sparse Unconstrained Optimisation Using Automatic Differentiation, Journal of
#' Optimization Theory and Applications, vol. 60, no. 2, pp. 261-275, 1989.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makeDixonPriceFunction = function(dimensions) {
  assertCount(dimensions)
  i = 1:dimensions
  force(dimensions)
  global.opt.params = 2^((-1) * (2^i - 2) / 2^i)
  makeSingleObjectiveFunction(
    name = paste(dimensions, "-d Dixon-Price function", sep = ""),
    id = paste0("dixonPrice_", dimensions, "d"),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      a = (x[1] - 1)^2
      i = 2:length(x)
      b = sum(i * (2 * x[i]^2 - x[i - 1])^2)
      return(a + b)
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(-10, dimensions),
      upper = rep(10, dimensions),
      vector = TRUE
    ),
    tags = attr(makeDixonPriceFunction, "tags"),
    global.opt.params = global.opt.params,
    global.opt.value = 0
  )
}

class(makeDixonPriceFunction) = c("function", "smoof_generator")
attr(makeDixonPriceFunction, "name") = c("Dixon-Price")
attr(makeDixonPriceFunction, "type") = c("single-objective")
attr(makeDixonPriceFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "scalable", "unimodal")
