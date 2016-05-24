#' Hansen Function
#'
#' Test function with multiple global optima based on the definition
#' \deqn{f(\mathbf{x}) = \sum_{i = 1}^{4} (i + 1)\cos(i\mathbf{x}_1 + i - 1) \sum_{j = 1}^{4} (j + 1)\cos((j + 2) \mathbf{x}_2 + j + 1)}
#' subject to \eqn{\mathbf{x}_i \in [-10, 10], i = 1, 2}.
#'
#' @references C. Fraley, Software Performances on Nonlinear lems, Technical
#' Report no. STAN-CS-89-1244, Computer Science, Stanford University, 1989.
#'
#' @template ret_smoof_single
#' @export
makeHansenFunction = function() {
  makeSingleObjectiveFunction(
    name = "Hansen Function",
    id = "hansen_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      i = j = 0:4
      a = sum((i + 1) * cos(i * x[1] + i + 1))
      b = sum((j + 1) * cos((j + 2) * x[2] + j + 1))
      return (a * b)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makeHansenFunction, "tags"),
    global.opt.params = matrix(
      c(-7.589893, -7.708314,
        -7.589893, -1.425128,
        -7.589893, 4.858057,
        -1.306708, -7.708314,
        -1.306708, 4.858057,
        4.976478, 4.858057,
        4.976478, -1.425128,
        4.976478, -7.708314),
      ncol = 2L, byrow = TRUE),
    global.opt.value = -176.5418
  )
}

class(makeHansenFunction) = c("function", "smoof_generator")
attr(makeHansenFunction, "name") = c("Hansen")
attr(makeHansenFunction, "type") = c("single-objective")
attr(makeHansenFunction, "tags") = c("single-objective", "continuous", "differentiable", "separable", "non-scalable", "multimodal")
