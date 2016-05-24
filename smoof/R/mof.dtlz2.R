#' DTLZ2 Function (family)
#'
#' Builds and returns the multi-objective DTLZ2 test problem.
#'
#' The DTLZ2 test problem is defined as follows:
#'
#' Minimize \eqn{f_1(\mathbf{x}) = (1+g(\mathbf{x}_M)) \cos(x_1\pi/2) \cos(x_2\pi/2) \cdots \cos(x_{M-2}\pi/2) \cos(x_{M-1}\pi/2),}{
#' f[1](X) = (1 + g(XM)) * cos(x[1] * pi/2) * cos(x[2] * pi/2) * ... * cos(x[M-2] * pi/2) * cos(x[M-1] * pi/2)}
#'
#' Minimize \eqn{f_2(\mathbf{x}) = (1+g(\mathbf{x}_M)) \cos(x_1\pi/2) \cos(x_2\pi/2) \cdots \cos(x_{M-2}\pi/2) \sin(x_{M-1}\pi/2),}{
#' f[2](X) = (1 + g(XM)) * cos(x[1] * pi/2) * cos(x[2] * pi/2) * ... * cos(x[M-2] * pi/2) * sin(x[M-1] * pi/2)}
#'
#' Minimize \eqn{f_3(\mathbf{x}) = (1+g(\mathbf{x}_M)) \cos(x_1\pi/2) \cos(x_2\pi/2) \cdots \sin(x_{M-2}\pi/2),}{
#' f[3](X) = (1 + g(XM)) * cos(x[1] * pi/2) * cos(x[2] * pi/2) * ... * sin(x[M-2] * pi/2)}
#'
#' \eqn{\vdots\\}{...}
#'
#' Minimize \eqn{f_{M-1}(\mathbf{x}) = (1+g(\mathbf{x}_M)) \cos(x_1\pi/2) \sin(x_2\pi/2),}{
#' f[M-1](X) = (1 + g(XM)) * cos(x[1] * pi/2) * sin(x[2] * pi/2)}
#'
#' Minimize \eqn{f_{M}(\mathbf{x}) = (1+g(\mathbf{x}_M)) \sin(x_1\pi/2),}{
#' f[M](X) = (1 + g(XM)) * sin(x[1] * pi/2)}
#'
#' with \eqn{0 \leq x_i \leq 1}{0 <= x[i] <= 1}, for \eqn{i=1,2,\dots,n,}{i=1,2,...,n}
#'
#' where \eqn{g(\mathbf{x}_M) = \sum\limits_{x_i\in \mathbf{x}_M}(x_i-0.5)^2}{
#' g(XM) = sum{x[i] in XM} {(x[i] - 0.5)^2}}
#'
#' @references K. Deb and L. Thiele and M. Laumanns and E. Zitzler. Scalable
#' Multi-Objective Optimization Test Problems. Computer Engineering and Networks
#' Laboratory (TIK), Swiss Federal Institute of Technology (ETH) Zurich, 112, 2001
#'
#' @param dimensions [\code{integer(1)}]\cr
#'   Number of decision variables.
#' @param n.objectives [\code{integer(1)}]\cr
#'   Number of objectives.
#' @return [\code{smoof_multi_objective_function}]
#' @export
makeDTLZ2Function = function(dimensions, n.objectives) {
  assertInt(n.objectives, lower = 2L, na.ok = FALSE)
  assertInt(dimensions, lower = n.objectives, na.ok = FALSE)

  # Renaming n.objectives here to stick to the notation in the paper
  M = n.objectives

  force(M)

  # C++ implementation
  fn = function(x) {
    assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
    dtlz_2(x, M)
  }

  # Reference R implementation
  # fn = function(x) {
  #   assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
  #   f = numeric(M)
  #   n = length(x)
  #   xm = x[M:n]
  #   g = sum((xm - 0.5)^2)
  #   a = (1 + g)
  #   prod.xi = 1
  #   for(i in M:2) {
  #     f[i] = a * prod.xi * sin(x[M - i + 1] * pi * 0.5)
  #     prod.xi = prod.xi * cos(x[M - i + 1] * pi * 0.5)
  #   }
  #   f[1] = a * prod.xi
  #   return(f)
  # }

  makeMultiObjectiveFunction(
    name = "DTLZ2 Function",
    id = paste0("dtlz2_", dimensions, "d_", n.objectives, "o"),
    description = "Deb et al.",
    fn = fn,
    par.set =  makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(0, dimensions),
      upper = rep(1, dimensions),
      vector = TRUE
      ),
    n.objectives = n.objectives,
    ref.point = rep(11, n.objectives)
  )
}

class(makeDTLZ2Function) = c("function", "smoof_generator")
attr(makeDTLZ2Function, "name") = c("DTLZ2")
attr(makeDTLZ2Function, "type") = c("multi-objective")
attr(makeDTLZ2Function, "tags") = c("multi-objective")
