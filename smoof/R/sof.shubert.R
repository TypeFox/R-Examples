#' Shubert Function
#'
#' The defintion of this two-dimensional function is given by
#' \deqn{f(\mathbf{x}) = \prod_{i = 1}^{2} \left(\sum_{j = 1}^{5} \cos((j + 1)\mathbf{x}_i + j\right)}
#' subject to \eqn{\mathbf{x}_i \in [-10, 10], i = 1, 2}.
#'
#' @references J. P. Hennart (ed.), Numerical Analysis, Proc. 3rd AS Workshop,
#' Lecture Notes in Mathematics, vol. 90, Springer, 1982.
#'
#' @template ret_smoof_single
#' @export
makeShubertFunction = function() {
  makeSingleObjectiveFunction(
    name = "Shubert Function",
    id = "shubert_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      j = 1:5
      a = sum(j * cos((j + 1) * x[1] + j))
      b = sum(j * cos((j + 1) * x[2] + j))
      return (a * b)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-10, -10),
      upper = c(10, 10),
      vector = TRUE
    ),
    tags = attr(makeShubertFunction, "tags"),
    global.opt.params = matrix(
      c(-7.0835, 4.8580,
        -7.0835, -7.7083,
        -1.4251, -7.0835,
        -1.4251, -0.8003,
        -7.7083, -7.0835,
        -7.7083, -0.8003,
        -0.8003, -7.7083,
        -0.8003, 4.8580,
        5.4828, -7.7083,
        5.4828, 4.8580,
        4.8580, 5.4828,
        -7.0835, -1.4251,
        -7.7083, 5.4828,
        -0.8003, -1.4251,
        -1.4251, 5.4828,
        4.8580,-7.0835,
        4.8580,-0.8003),
      ncol = 2L, byrow = TRUE),
    global.opt.value = -186.7309
  )
}

class(makeShubertFunction) = c("function", "smoof_generator")
attr(makeShubertFunction, "name") = c("Shubert")
attr(makeShubertFunction, "type") = c("single-objective")
attr(makeShubertFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-scalable", "multimodal")
