#' Bent-Cigar Function
#'
#' Scalable test function \eqn{f} with
#' \deqn{f(\mathbf{x}) = x_1^2 + 10^6 \sum_{i = 2}^{n} x_i^2}
#' subject to \eqn{-100 \leq \mathbf{x}_i \leq 100} for \eqn{i = 1, \ldots, n}.
#'
#' @references See \url{http://al-roomi.org/benchmarks/unconstrained/n-dimensions/164-bent-cigar-function}.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makeBentCigarFunction = function(dimensions) {
  assertCount(dimensions, na.ok = FALSE)
  force(dimensions)
  makeSingleObjectiveFunction(
    name = "Bent-Cigar Function",
    id = "bentCigar_2d",
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      x[1]^2 + 1e+06 * sum(x[2:dimensions]^2)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-100, -100),
      upper = c(100, 100),
      vector = TRUE
    ),
    tags = attr(makeBentCigarFunction, "tags"),
    global.opt.params = rep(0, dimensions),
    global.opt.value = 0
  )
}

class(makeBentCigarFunction) = c("function", "smoof_generator")
attr(makeBentCigarFunction, "name") = c("Bent-Cigar")
attr(makeBentCigarFunction, "type") = c("single-objective")
attr(makeBentCigarFunction, "tags") = c("single-objective", "continuous", "differentiable", "non-separable", "scalable", "unimodal")
