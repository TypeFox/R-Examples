#' Alpine02 function
#'
#' Another multimodal optimization test function. The implementation is based on
#' the formula
#' \deqn{f(\mathbf{x}) = \prod_{i = 1}^{n} \sqrt{\mathbf{x}_i}\sin(\mathbf{x}_i)}
#' with \eqn{\mathbf{x}_i \in [0, 10]} for \eqn{i = 1, \ldots, n}.
#'
#' @references M. Clerc, The Swarm and the Queen, Towards a Deterministic and
#' Adaptive Particle Swarm Optimization, IEEE Congress on Evolutionary Computation,
#' Washington DC, USA, pp. 1951-1957, 1999.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makeAlpine02Function = function(dimensions) {
  assertCount(dimensions)
  force(dimensions)
  makeSingleObjectiveFunction(
    name = paste(dimensions, "-d Alpine N. 2 Function", sep = ""),
    id = paste0("alpine02_", dimensions, "d"),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      prod(sqrt(x) * sin(x))
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(0, dimensions),
      upper = rep(10, dimensions),
      vector = TRUE
    ),
    minimize = FALSE,
    tags = attr(makeAlpine02Function, "tags"),
    global.opt.params = rep(7.917, dimensions),
    global.opt.value = 2.808^dimensions
  )
}

class(makeAlpine02Function) = c("function", "smoof_generator")
attr(makeAlpine02Function, "name") = c("Alpine N. 2")
attr(makeAlpine02Function, "type") = c("single-objective")
attr(makeAlpine02Function, "tags") = c("single-objective", "continuous", "differentiable", "separable", "scalable", "multimodal")
