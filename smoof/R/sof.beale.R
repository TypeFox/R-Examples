#' Beale Function
#'
#' Multimodal single-objective test function for optimization. It is based on
#' the mathematic formula
#' \deqn{f(\mathbf{x}) = (1.5 - \mathbf{x}_1 + \mathbf{x}_1\mathbf{x}_2)^2 + (2.25 - \mathbf{x}_1 + \mathbf{x}_1\mathbf{x}_2^2)^2 + (2.625 - \mathbf{x}_1 + \mathbf{x}_1\mathbf{x}_2^3)^2}
#' usually evaluated within the bounds \eqn{\mathbf{x}_i \in [-4.5, 4.5], i = 1, 2}.
#' The function has a flat but multimodal region aroung the single global optimum
#' and large peaks in the edges of its definition space.
#'
#' @template ret_smoof_single
#' @export
makeBealeFunction = function() {
  makeSingleObjectiveFunction(
    name = "Beale Function",
    id = "beale_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      a = x[1] * x[2]
      b = a * x[2]
      c = b * x[2]
      (1.5 - x[1] + a)^2 + (2.25 - x[1] + b)^2 + (2.625 - x[1] + c)^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-4.5, -4.5),
      upper = c(4.5, 4.5),
      vector = TRUE
    ),
    tags = attr(makeBealeFunction, "tags"),
    global.opt.params = c(3, 0.5),
    global.opt.value = 0
  )
}

class(makeBealeFunction) = c("function", "smoof_generator")
attr(makeBealeFunction, "name") = c("Beale")
attr(makeBealeFunction, "type") = c("single-objective")
attr(makeBealeFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "unimodal")
