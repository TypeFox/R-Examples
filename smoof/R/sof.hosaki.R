#' Hosaki Function
#'
#' Two-dimensional test function \eqn{f} with
#' \deqn{f(\mathbf{x}) = (1 - 8 \mathbf{x}_1 + 7 \mathbf{x}_1^2 - 7/3 \mathbf{x}_1^3 + 1/4 \mathbf{x}_1^4)\mathbf{x}_2^2e^{-\mathbf{x}_2}}
#' subject to \eqn{0 \leq \mathbf{x}_1 \leq 5} and \eqn{0 \leq \mathbf{x}_2 \leq 6}.
#'
#' @references G. A. Bekey, M. T. Ung, A Comparative Evaluation of Two Global
#' Search Algorithms, IEEE Transaction on Systems, Man and Cybernetics, vol. 4,
#' no. 1, pp. 112- 116, 1974.
#'
#' @template ret_smoof_single
#' @export
makeHosakiFunction = function() {
  makeSingleObjectiveFunction(
    name = "Hosaki Function",
    id = "hosaki_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      (1 - 8 * x[1] + 7 * x[1]^2 - 7 * x[1]^3 / 3 + 0.25 * x[1]^4) * x[2]^2 * exp(-x[2])
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(0, 0),
      upper = c(5, 6),
      vector = TRUE
    ),
    tags = attr(makeHosakiFunction, "tags"),
    global.opt.params = c(4, 2),
    global.opt.value = -2.3458
  )
}

class(makeHosakiFunction) = c("function", "smoof_generator")
attr(makeHosakiFunction, "name") = c("Hosaki")
attr(makeHosakiFunction, "type") = c("single-objective")
attr(makeHosakiFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "multimodal")
