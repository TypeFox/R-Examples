#' DTLZ5 Function (family)
#'
#' Builds and returns the multi-objective DTLZ5 test problem. This problem
#' can be characterized by a disconnected Pareto-optimal front in the search
#' space. This introduces a new challenge to evolutionary multi-objective
#' optimizers, i.e., to maintain different subpopulations within the search
#' space to cover the entire Pareto-optimal front.
#'
#' The DTLZ5 test problem is defined as follows:
#'
#' Minimize \eqn{f_1(\mathbf{x}) = (1+g(\mathbf{x}_M)) \cos(\theta_1\pi/2) \cos(\theta_2\pi/2) \cdots \cos(\theta_{M-2}\pi/2) \cos(\theta_{M-1}\pi/2),}{
#' f[1](X) = (1 + g(XM)) * cos(theta[1] * pi/2) * cos(theta[2] * pi/2) * ... * cos(theta[M-2] * pi/2) * cos(theta[M-1] * pi/2)}
#'
#' Minimize \eqn{f_2(\mathbf{x}) = (1+g(\mathbf{x}_M)) \cos(\theta_1\pi/2) \cos(\theta_2\pi/2) \cdots \cos(\theta_{M-2}\pi/2) \sin(\theta_{M-1}\pi/2),}{
#' f[2](X) = (1 + g(XM)) * cos(theta[1] * pi/2) * cos(theta[2] * pi/2) * ... * cos(theta[M-2] * pi/2) * sin(theta[M-1] * pi/2)}
#'
#' Minimize \eqn{f_3(\mathbf{x}) = (1+g(\mathbf{x}_M)) \cos(\theta_1\pi/2) \cos(\theta_2\pi/2) \cdots \sin(\theta_{M-2}\pi/2),}{
#' f[3](X) = (1 + g(XM)) * cos(theta[1] * pi/2) * cos(theta[2] * pi/2) * ... * sin(theta[M-2] * pi/2)}
#'
#' \eqn{\vdots\\}{...}
#'
#' Minimize \eqn{f_{M-1}(\mathbf{x}) = (1+g(\mathbf{x}_M)) \cos(\theta_1\pi/2) \sin(\theta_2\pi/2),}{
#' f[M-1](X) = (1 + g(XM)) * cos(theta[1] * pi/2) * sin(theta[2] * pi/2)}
#'
#' Minimize \eqn{f_{M}((1+g(\mathbf{x}_M)) \sin(\theta_1\pi/2),}{
#' f[M](X) = (1 + g(XM)) * sin(theta[1] * pi/2)}
#'
#' with \eqn{0 \leq x_i \leq 1}{0 <= x[i] <= 1}, for \eqn{i=1,2,\dots,n,}{i=1,2,...,n}
#'
#' where \eqn{\theta_i = \frac{\pi}{4(1+ g(\mathbf{x}_M))} (1+2g(\mathbf{x}_M)x_i),}{
#' theta[i] = pi / (4 * (1 + g(XM))) * (1 + 2 * g(XM) * x[i]),}
#' for \eqn{i = 2,3,\dots,(M-1)}{i = 2,3,...,(M-1)}
#'
#' and \eqn{g(\mathbf{x}_M) = \sum\limits_{x_i\in\mathbf{x}_M}(x_i-0.5)^2}{
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
makeDTLZ5Function = function(dimensions, n.objectives) {
  assertInt(n.objectives, lower = 2L, na.ok = FALSE)
  assertInt(dimensions, lower = n.objectives, na.ok = FALSE)

  # Renaming n.objectives here to stick to the notation in the paper
  M = n.objectives

  force(M)

  # C++ implementation
  fn = function(x) {
    assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
    dtlz_5(x, M)
  }

  # Reference R implementation
  # fn = function(x) {
  #   assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
  #   f = numeric(M)
  #   n = length(x)
  #   theta = numeric(M-1)
  #   xm = x[M:n]
  #   g = sum((xm - 0.5)^2)
  #   t = pi / (4 * (1 + g))
  #   theta[1] = x[1] * pi / 2
  #   theta[-1] = t * (1 + 2 * g * x[2:(M - 1)])
  #   a = (1 + g)
  #   prod.xi = 1
  #   for(i in M:2) {
  #     f[i] = a * prod.xi * sin(theta[M - i + 1])
  #     prod.xi = prod.xi * cos(theta[M - i + 1])
  #   }
  #   f[1] = a * prod.xi
  #   return(f)
  # }

  makeMultiObjectiveFunction(
    name = "DTLZ5 Function",
    id = paste0("dtlz5_", dimensions, "d_", n.objectives, "o"),
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

class(makeDTLZ5Function) = c("function", "smoof_generator")
attr(makeDTLZ5Function, "name") = c("DTLZ5")
attr(makeDTLZ5Function, "type") = c("multi-objective")
attr(makeDTLZ5Function, "tags") = c("multi-objective")
