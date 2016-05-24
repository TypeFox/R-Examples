#' Styblinkski-Tang function
#'
#' This function is based on the defintion
#' \deqn{f(\mathbf{x}) = \frac{1}{2} \sum_{i = 1}^{2} (\mathbf{x}_i^4 - 16 \mathbf{x}_i^2 + 5\mathbf{x}_i)}
#' with box-constraints given by \eqn{\mathbf{x}_i \in [-5, 5], i = 1, 2}.
#'
#' @references Z. K. Silagadze, Finding Two-Dimesnional Peaks, Physics of
#' Particles and Nuclei Letters, vol. 4, no. 1, pp. 73-80, 2007.
#'
#' @template ret_smoof_single
#' @export
#FIXME: hmm, this can be formulated as a scalable problem
makeStyblinkskiTangFunction = function() {
  makeSingleObjectiveFunction(
    name = "Styblinkski-Tang Function",
    id = "styblinskiTang_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      a = x^2
      b = a^2
      return(0.5 * sum(b - 16 * a + 5 * x))
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = rep(-5, 2L),
      upper = rep(5, 2L),
      vector = TRUE
    ),
    tags = attr(makeStyblinkskiTangFunction, "tags"),
    global.opt.params = c(-2.903534, -2.903534),
    global.opt.value = -78.332
  )
}

class(makeStyblinkskiTangFunction) = c("function", "smoof_generator")
attr(makeStyblinkskiTangFunction, "name") = c("Styblinkski-Tang")
attr(makeStyblinkskiTangFunction, "type") = c("single-objective")
attr(makeStyblinkskiTangFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
