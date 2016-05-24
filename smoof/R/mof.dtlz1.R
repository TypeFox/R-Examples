#' DTLZ1 Function (family)
#'
#' Builds and returns the multi-objective DTLZ1 test problem.
#'
#' The DTLZ1 test problem is defined as follows:
#'
#' Minimize \eqn{f_1(\mathbf{x}) = \frac{1}{2} x_1 x_2 \cdots x_{M-1} (1+g(\mathbf{x}_M),}{
#' f[1](X) = 1/2 * x[1] * x[2] * ... * x[M-1] * (1 + g(XM))}
#'
#' Minimize \eqn{f_2(\mathbf{x}) = \frac{1}{2} x_1 x_2 \cdots (1-x_{M-1}) (1+g(\mathbf{x}_M)),}{
#' f[2](X) = 1/2 * x[1] * x[2] * ... * (1 - x[M-1]) * (1 + g(XM))}
#'
#' \eqn{\vdots\\}{...}
#'
#' Minimize \eqn{f_{M-1}(\mathbf{x}) = \frac{1}{2} x_1 (1-x_2) (1+g(\mathbf{x}_M)),}{
#' f[M-1](X) = 1/2 * x[1] * (1 - x[2]) * (1 + g(XM))}
#'
#' Minimize \eqn{f_{M}(\mathbf{x}) = \frac{1}{2} (1-x_1) (1+g(\mathbf{x}_M)),}{
#' f[M](X) = 1/2 * (1 - x[1]) * (1 + g(XM))}
#'
#' with \eqn{0 \leq x_i \leq 1}{0 <= x[i] <= 1}, for \eqn{i=1,2,\dots,n,}{i=1,2,...,n}
#'
#' where \eqn{g(\mathbf{x}_M) = 100 \left[|\mathbf{x}_M| + \sum\limits_{x_i \in \mathbf{x}_M} (x_i - 0.5)^2 - \cos(20\pi(x_i - 0.5))\right]}{
#' g(XM) = 100 * (|XM| + sum{x[i] in XM} {(x[i] - 0.5)^2 - cos(20 * pi * (x[i] - 0.5))})}
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
makeDTLZ1Function = function(dimensions, n.objectives) {
  assertInt(n.objectives, lower = 2L, na.ok = FALSE)
  assertInt(dimensions, lower = n.objectives, na.ok = FALSE)

  # Renaming vars here to stick to the notation in the paper
  # number of decision variables in the last group (see x_m in the paper)
  k = dimensions - n.objectives + 1
  M = n.objectives

  force(M)
  force(k)

  # C++ implementation
  fn = function(x) {
    assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
    dtlz_1(x, M)
  }

  # Reference R implementation
  # fn = function(x) {
  #   assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
  #   f = numeric(M)
  #   n = length(x)
  #   xm = x[M:n]
  #   g = 100 * (k + sum((xm - 0.5)^2 - cos(20 * pi * (xm - 0.5))))
  #   a = 0.5 * (1 + g)
  #   prod.xi = 1
  #   for(i in M:2) {
  #     f[i] = a * prod.xi * (1 - x[M - i + 1])
  #     prod.xi = prod.xi * x[M - i + 1]
  #   }
  #   f[1] = a * prod.xi
  #   return(f)
  # }

  makeMultiObjectiveFunction(
    name = "DTLZ1 Function",
    id = paste0("dtlz1_", dimensions, "d_", n.objectives, "o"),
    description = "Deb et al.",
    fn = fn,
    par.set =  makeNumericParamSet(
      len = dimensions,
      id = "x",
      #FIXME: any box constraints?
      lower = rep(0, dimensions),
      upper = rep(1, dimensions),
      vector = TRUE
      ),
    n.objectives = n.objectives,
    ref.point = rep(11, n.objectives)
  )
}

class(makeDTLZ1Function) = c("function", "smoof_generator")
attr(makeDTLZ1Function, "name") = c("DTLZ1")
attr(makeDTLZ1Function, "type") = c("multi-objective")
attr(makeDTLZ1Function, "tags") = c("multi-objective")
