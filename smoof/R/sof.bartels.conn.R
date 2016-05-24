#' Bartels Conn Function
#'
#' The Bartels Conn Function is defined as
#' \deqn{f(\mathbf{x}) = |\mathbf{x}_1^2 + \mathbf{x}_2^2 + \mathbf{x}_1\mathbf{x}_2| + |\sin(\mathbf{x}_1)| + |\cos(\mathbf{x})|}
#' subject to \eqn{\mathbf{x}_i \in [-500, 500]} for \eqn{i = 1, 2}.
#'
#' @template ret_smoof_single
#' @export
makeBartelsConnFunction = function() {
  makeSingleObjectiveFunction(
    name = "Bartels Conn Function",
    id = "bartelsConn_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      abs(x[1]^2 + x[2]^2 + x[1] * x[2]) + abs(sin(x[1])) + abs(cos(x[2]))
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-500, -500),
      upper = c(500, 500),
      vector = TRUE
    ),
    tags = attr(makeBartelsConnFunction, "tags"),
    global.opt.params = c(0, 0),
    global.opt.value = 1
  )
}

class(makeBartelsConnFunction) = c("function", "smoof_generator")
attr(makeBartelsConnFunction, "name") = c("Bartels Conn")
attr(makeBartelsConnFunction, "type") = c("single-objective")
attr(makeBartelsConnFunction, "tags") = c("single-objective", "continuous", "non-differentiable", "non-separable", "non-scalable", "multimodal")
