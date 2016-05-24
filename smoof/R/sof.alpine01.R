#' Alpine01 function
#'
#' Highly multimodal single-objective optimization test function. It is defined
#' as \deqn{f(\mathbf{x}) = \sum_{i = 1}^{n} |\mathbf{x}_i \sin(\mathbf{x}_i) + 0.1\mathbf{x}_i|}
#' with box constraints \eqn{\mathbf{x}_i \in [-10, 10]} for \eqn{i = 1, \ldots, n}.
#'
#' @references S. Rahnamyan, H. R. Tizhoosh, N. M. M. Salama, A Novel Population
#' Initialization Method for Accelerating Evolutionary Algorithms, Computers and
#' Mathematics with Applications, vol. 53, no. 10, pp. 1605-1614, 2007.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makeAlpine01Function = function(dimensions) {
  assertCount(dimensions)
  force(dimensions)
  makeSingleObjectiveFunction(
    name = paste(dimensions, "-d Alpine01 Function", sep = ""),
    id = paste0("alpine01_", dimensions, "d"),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      sum(abs(x * sin(x) + 0.1 * x))
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(-10, dimensions),
      upper = rep(10, dimensions),
      vector = TRUE
    ),
    tags = attr(makeAlpine01Function, "tags"),
    global.opt.params = rep(0, dimensions),
    global.opt.value = 0
  )
}

class(makeAlpine01Function) = c("function", "smoof_generator")
attr(makeAlpine01Function, "name") = c("Alpine N. 1")
attr(makeAlpine01Function, "type") = c("single-objective")
attr(makeAlpine01Function, "tags") = c("single-objective", "continuous", "non-differentiable", "separable", "scalable", "multimodal")
