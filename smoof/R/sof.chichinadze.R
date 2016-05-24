#' Chichinadze Function
#'
#' Continuous single-objective test function \eqn{f} with
#' \deqn{f(\mathbf{x}) = \mathbf{x}_1^2 - 12 \mathbf{x}_1 + 11 + 10\cos(0.5\pi\mathbf{x}_1) + 8\sin(2.5\pi\mathbf{x}_1) - (0.25)^{0.5}\exp(-0.5(\mathbf{x}_2 - 0.5)^2)}
#' with \eqn{-30 \leq \mathbf{x}_i \leq 30, i = 1, 2}.
#'
#' @template ret_smoof_single
#' @export
makeChichinadzeFunction = function() {
  makeSingleObjectiveFunction(
    name = "Chichinadze Function",
    id = "chichinadze_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      x[1]^2 - 12 * x[1] + 11 + 10 * cos(pi * 0.5 * x[1]) + 8 * sin(5 * pi * x[1]) - exp(-0.5 * (x[2] - 0.5)^2) / sqrt(5)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-30, -30),
      upper = c(30, 30),
      vector = TRUE
    ),
    tags = attr(makeChichinadzeFunction, "tags"),
    global.opt.params = c(5.90133, 0.5),
    global.opt.value = -43.3159
  )
}

class(makeChichinadzeFunction) = c("function", "smoof_generator")
attr(makeChichinadzeFunction, "name") = c("Chichinadze")
attr(makeChichinadzeFunction, "type") = c("single-objective")
attr(makeChichinadzeFunction, "tags") = c("single-objective", "continuous", "differentiable", "separable", "non-scalable", "multimodal")
