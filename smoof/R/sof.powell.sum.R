#' Powell-Sum Function
#'
#' The formula that underlies the implementation is given by
#' \deqn{f(\mathbf{x}) = \sum_{i=1}^n |\mathbf{x}_i|^{i+1}}
#' with \eqn{\mathbf{x}_i \in [-1, 1], i = 1, \ldots, n}.
#'
#' @references S. Rahnamyan, H. R. Tizhoosh, N. M. M. Salama, A Novel Population
#' Initialization Method for Accelerating Evolutionary Algorithms, Computers and
#' Mathematics with Applications, vol. 53, no. 10, pp. 1605-1614, 2007.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makePowellSumFunction = function(dimensions) {
  assertCount(dimensions)
  force(dimensions)
  makeSingleObjectiveFunction(
    name = paste(dimensions, "-d Powell-Sum Function", sep = ""),
    id = paste0("powellSum", dimensions, "d"),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      a = (1:length(x)) + 1L
      sum(abs(x)^a)
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(-1, dimensions),
      upper = rep(1, dimensions),
      vector = TRUE
    ),
    tags = attr(makePowellSumFunction, "tags"),
    global.opt.params = rep(0, dimensions),
    global.opt.value = 0
  )
}

class(makePowellSumFunction) = c("function", "smoof_generator")
attr(makePowellSumFunction, "name") = c("Double-Sum")
attr(makePowellSumFunction, "type") = c("single-objective")
attr(makePowellSumFunction, "tags") = c("single-objective", "continuous", "differentiable", "separable", "scalable", "unimodal")
