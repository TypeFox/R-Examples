#' @title Judge function.
#'
#' @description Two-dimensional test function based on the formula
#' \deqn{f(\mathbf{x}) = \sum_{i=1}^{20} \left[(x_1 + B_i x_2 + C_i x_2^2) - A_i\right]^2}
#' with \eqn{\mathbf{x}_1, \mathbf{x}_2 \in [-10, 10]}. For details on \eqn{A, B, C}
#' see the referenced website.
#'
#' @references See \url{http://al-roomi.org/benchmarks/unconstrained/2-dimensions/133-judge-s-function}.
#'
#' @template ret_smoof_single
#' @export
makeJudgeFunction = function() {
  makeSingleObjectiveFunction(
    name = "Judge Function",
    id = "judge_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      A = c(4.284, 4.149, 3.877, 0.533, 2.211, 2.389, 2.145, 3.231, 1.998, 1.379, 2.106, 1.428, 1.011, 2.179, 2.858, 1.388, 1.651, 1.593, 1.046, 2.152)
      B = c(0.286, 0.973, 0.384, 0.276, 0.973, 0.543, 0.957, 0.948, 0.543, 0.797, 0.936, 0.889, 0.006, 0.828, 0.399, 0.617, 0.939, 0.784, 0.072, 0.889)
      C = c(0.645, 0.585, 0.310, 0.058, 0.455, 0.779, 0.259, 0.202, 0.028, 0.099, 0.142, 0.296, 0.175, 0.180, 0.842, 0.039, 0.103, 0.620, 0.158, 0.704)
      i = 1:20
      sum(((x[1] + B * x[2] + C * x[2]^2) - A)^2)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makeJudgeFunction, "tags"),
    global.opt.params = c(0.864787285816574, 1.235748499036571),
    global.opt.value = 16.081730132960381
  )
}

class(makeJudgeFunction) = c("function", "smoof_generator")
attr(makeJudgeFunction, "name") = c("Judge")
attr(makeJudgeFunction, "type") = c("single-objective")
attr(makeJudgeFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
