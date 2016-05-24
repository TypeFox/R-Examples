#' Cross-In-Tray Function
#'
#' Non-scalable, two-dimensional test function for numerical optimization with
#' \deqn{f(\mathbf{x}) = -0.0001\left(|\sin(\mathbf{x}_1\mathbf{x}_2\exp(|100 - [(\mathbf{x}_1^2 + \mathbf{x}_2^2)]^{0.5} / \pi|)| + 1\right)^{0.1}}
#' subject to \eqn{\mathbf{x}_i \in [-15, 15]} for \eqn{i = 1, 2}.
#'
#' @references S. K. Mishra, Global Optimization By Differential Evolution and
#' Particle Swarm Methods: Evaluation On Some Benchmark Functions, Munich
#' Research Papers in Economics.
#'
#' @template ret_smoof_single
#' @export
makeCrossInTrayFunction = function() {
  makeSingleObjectiveFunction(
    name = "Cross-In-Tray Function",
    id = "crossInTray_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      a = exp(abs(100 - (sqrt(x[1]^2 + x[2]^2) / pi)))
      -0.0001 * (abs(a * sin(x[1]) * sin(x[2])) + 1)^(0.1)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-15, -15),
      upper = c(15, 15),
      vector = TRUE
    ),
    tags = attr(makeCrossInTrayFunction, "tags"),
    global.opt.params = matrix(
      c(1.349406685353340, 1.349406608602084,
        1.349406685353340, -1.349406608602084,
        -1.349406685353340, 1.349406608602084,
        -1.349406685353340, -1.349406608602084),
      ncol = 2L, byrow = TRUE),
    global.opt.value = -2.06261218
  )
}

class(makeCrossInTrayFunction) = c("function", "smoof_generator")
attr(makeCrossInTrayFunction, "name") = c("Cross-In-Tray")
attr(makeCrossInTrayFunction, "type") = c("single-objective")
attr(makeCrossInTrayFunction, "tags") = c("single-objective", "continuous", "non-separable", "non-scalable", "multimodal")
