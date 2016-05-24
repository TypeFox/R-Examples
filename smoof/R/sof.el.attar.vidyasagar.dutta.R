#' El-Attar-Vidyasagar-Dutta Function
#'
#' This function is based on the formula
#' \deqn{f(\mathbf{x}) = (\mathbf{x}_1^2 + \mathbf{x}_2 - 10)^2 + (\mathbf{x}_1 + \mathbf{x}_2^2 - 7)^2 + (\mathbf{x}_1^2 + \mathbf{x}_2^3 - 1)^2}
#' subject to \eqn{\mathbf{x}_i \in [-500, 500], i = 1, 2}.
#'
#' @references R. A. El-Attar, M. Vidyasagar, S. R. K. Dutta, An Algorithm for
#' II-norm Minimiza- tion With Application to Nonlinear II-approximation, SIAM
#' Journal on Numverical Analysis, vol. 16, no. 1, pp. 70-86, 1979.
#'
#' @template ret_smoof_single
#' @export
makeElAttarVidyasagarDuttaFunction = function() {
  makeSingleObjectiveFunction(
    name = "El-Attar-Vidyasagar-Dutta Function",
    id = "elAttarVidyasagarDutta_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      (x[1]^2 + x[2] - 10)^2 + (x[1] + x[2]^2 - 7)^2 + (x[1]^2 + x[2]^3 - 1)^2
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-100, -100),
      upper = c(100, 100),
      vector = TRUE
    ),
    tags = attr(makeElAttarVidyasagarDuttaFunction, "tags"),
    global.opt.params = c(3.40918683, -2.17143304),
    global.opt.value = 1.712780354
  )
}

class(makeElAttarVidyasagarDuttaFunction) = c("function", "smoof_generator")
attr(makeElAttarVidyasagarDuttaFunction, "name") = c("El-Attar-Vidyasagar-Dutta")
attr(makeElAttarVidyasagarDuttaFunction, "type") = c("single-objective")
attr(makeElAttarVidyasagarDuttaFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "non-scalable", "unimodal")
