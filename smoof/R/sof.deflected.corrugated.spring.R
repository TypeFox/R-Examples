#' @title Deflected Corrugated Spring function
#'
#' @description Scalable single-objective test function based on the formula
#' \deqn{f(\mathbf{x}) = 0.1 \sum_{i = 1}^{n} (x_i - \alpha)^2 - \cos\left(K \sqrt{\sum_{i = 1}^{n} (x_i - \alpha)^2}\right)}
#' with \eqn{\mathbf{x}_i \in [0, 2\alpha], i = 1, \ldots, n} and \eqn{\alpha = K = 5}
#' by default.
#'
#' @references See \url{http://al-roomi.org/benchmarks/unconstrained/n-dimensions/238-deflected-corrugated-spring-function}.
#'
#' @template arg_dimensions
#' @param K [\code{numeric(1)}]\cr
#'   Parameter.
#'   Default is 5.
#' @param alpha [\code{numeric(1)}]\cr
#'   Parameter.
#'   Default is 5.
#' @template ret_smoof_single
#' @export
makeDeflectedCorrugatedSpringFunction = function(dimensions, K = 5, alpha = 5) {
  assertCount(dimensions)
  assertNumber(K, na.ok = FALSE)
  assertNumber(alpha, na.ok = FALSE)

  force(K)
  force(alpha)
  force(dimensions)

  makeSingleObjectiveFunction(
    name = paste(dimensions, "-d Deflected Corrugated Spring function", sep = ""),
    #FIXME: eventually encode K and alpha in fun?
    id = paste0("deflectedCorrugatedSpring_", dimensions, "d_K", K, "alpha", alpha),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      a = (x - alpha)^2
      sa = sum(a)
      0.1 * sa - cos(K * sqrt(sa))
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(0, dimensions),
      upper = rep(2 * alpha, dimensions),
      vector = TRUE
    ),
    tags = attr(makeDeflectedCorrugatedSpringFunction, "tags"),
    global.opt.params = rep(alpha, dimensions),
    global.opt.value = -1
  )
}

class(makeDeflectedCorrugatedSpringFunction) = c("function", "smoof_generator")
attr(makeDeflectedCorrugatedSpringFunction, "name") = c("Deflected Corrugated Spring")
attr(makeDeflectedCorrugatedSpringFunction, "type") = c("single-objective")
attr(makeDeflectedCorrugatedSpringFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "scalable", "multimodal")
