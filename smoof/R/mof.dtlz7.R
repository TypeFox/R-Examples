#' DTLZ7 Function (family)
#'
#' Builds and returns the multi-objective DTLZ7 test problem. This problem
#' can be characterized by a disconnected Pareto-optimal front in the search
#' space. This introduces a new challenge to evolutionary multi-objective
#' optimizers, i.e., to maintain different subpopulations within the search
#' space to cover the entire Pareto-optimal front.
#'
#' The DTLZ7 test problem is defined as follows:
#'
#' Minimize \eqn{f_1(\mathbf{x}) = x_1,}{
#' f[1](X) = 1/2 * x[1] * x[2] * ... * x[M-1] * (1 + g(XM))}
#'
#' Minimize \eqn{f_2(\mathbf{x}) = x_2,}{
#' f[2](X) = 1/2 * x[1] * x[2] * ... * (1 - x[M-1]) * (1 + g(XM))}
#'
#' \eqn{\vdots\\}{...}
#'
#' Minimize \eqn{f_{M-1}(\mathbf{x}) = x_{M-1},}{
#' f[M-1](X) = 1/2 * x[1] * (1 - x[2]) * (1 + g(XM))}
#'
#' Minimize \eqn{f_{M}(\mathbf{x}) = (1+g(\mathbf{x}_M)) h(f_1,f_2,\cdots,f_{M-1}, g),}{
#' f[M](X) = 1/2 * (1 - x[1]) * (1 + g(XM))}
#'
#' with \eqn{0 \leq x_i \leq 1}{0 <= x[i] <= 1}, for \eqn{i=1,2,\dots,n,}{i=1,2,...,n}
#'
#' where \eqn{g(\mathbf{x}_M) = 1 + \frac{9}{|\mathbf{x}_M|} \sum_{x_i\in\mathbf{x}_M} x_i}{
#' g(XM) = 1 + 9 / |XM| * sum{x[i] in XM} {x[i]}}
#'
#' and \eqn{h(f_1,f_2,\cdots,f_{M-1}, g) = M - \sum_{i=1}^{M-1}\left[\frac{f_i}{1+g}(1 + sin(3\pi f_i))\right]}{
#' h(f[1],f[2],...f[M-1],g) = M - sum{i in 1:(M-1)} {f[i] / (1 + g) * (1 + sin(3 * pi * f[i]))}}
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
makeDTLZ7Function = function(dimensions, n.objectives) {
  assertInt(n.objectives, lower = 2L, na.ok = FALSE)
  assertInt(dimensions, lower = n.objectives, na.ok = FALSE)

  # Renaming n.objectives here to stick to the notation in the paper
  M = n.objectives

  force(M)

  # C++ implementation
  fn = function(x) {
    assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
    dtlz_7(x, M)
  }

  # Reference R implementation
  # fn = function(x) {
  #   f = numeric(M)
  #   n = length(x)
  #   f[1:(M - 1)] = x[1:(M - 1)]
  #   xm = x[(n - k):n]
  #   g = 1 + 9 * sum(xm) / k
  #   fi = f[1:(M - 1)]
  #   h = M - sum(fi  * (1 + sin(3 * pi * fi)) / (1 + g))
  #   f[M] = (1 + g) * h
  #   return(f)
  # }

  makeMultiObjectiveFunction(
    name = "DTLZ7 Function",
    id = paste0("dtlz7_", dimensions, "d_", n.objectives, "o"),
    description = "Deb et al.",
    fn = fn,
    par.set =  makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = rep(0, dimensions),
      upper = rep(1, dimensions),
      vector = TRUE
    ),
    n.objectives = n.objectives
  )
}

class(makeDTLZ7Function) = c("function", "smoof_generator")
attr(makeDTLZ7Function, "name") = c("DTLZ7")
attr(makeDTLZ7Function, "type") = c("multi-objective")
attr(makeDTLZ7Function, "tags") = c("multi-objective")
