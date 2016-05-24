#' Chung Reynolds Function
#'
#' The defintion is given by
#' \deqn{f(\mathbf{x}) = \left(\sum_{i=1}^{n} \mathbf{x}_i^2\right)^2}
#' with box-constraings \eqn{\mathbf{x}_i \in [-100, 100], i = 1, \ldots, n}.
#'
#' @references C. J. Chung, R. G. Reynolds, CAEP: An Evolution-Based Tool for
#' Real-Valued Func- tion Optimization Using Cultural Algorithms, International
#' Journal on Artificial In- telligence Tool, vol. 7, no. 3, pp. 239-291, 1998.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makeChungReynoldsFunction = function(dimensions) {
  assertCount(dimensions)
  force(dimensions)
  makeSingleObjectiveFunction(
    name = paste(dimensions, "-d Chung Reynolds Function", sep = ""),
    id = paste0("chungReynolds", dimensions, "d"),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      sum(x^2)^2
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(-100, dimensions),
      upper = rep(100, dimensions),
      vector = TRUE
    ),
    tags = attr(makeChungReynoldsFunction, "tags"),
    global.opt.params = rep(0, dimensions),
    global.opt.value = 0
  )
}

class(makeChungReynoldsFunction) = c("function", "smoof_generator")
attr(makeChungReynoldsFunction, "name") = c("Chung Reynolds")
attr(makeChungReynoldsFunction, "type") = c("single-objective")
attr(makeChungReynoldsFunction, "tags") = c("single-objective", "unimodal", "continuous", "differentiable", "scalable")
