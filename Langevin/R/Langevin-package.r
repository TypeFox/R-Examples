#' An \R package for stochastic data analysis
#'
#' The \pkg{Langevin} package provides functions to estimate drift and
#' diffusion functions from data sets.
#'
#'
#' This package was developed by the research group
#' \emph{Turbulence, Wind energy and Stochastics} (TWiSt) at the Carl von
#' Ossietzky University of Oldenburg (Germany).
#'
#' @section Mathematical Background: A wide range of dynamic systems can be
#' described by a stochastic differential equation, the Langevin equation. The
#' time derivative of the system trajectory \eqn{\dot{X}(t)} can be expressed as
#' a sum of a deterministic part \eqn{D^{(1)}} and the product of a stochastic
#' force \eqn{\Gamma(t)} and a weight coefficient \eqn{D^{(2)}}. The stochastic
#' force \eqn{\Gamma(t)} is \eqn{\delta}-correlated Gaussian white noise.
#'
#' For stationary continuous Markov processes Siegert et al. and Friedrich et
#' al. developed a method to reconstruct drift \eqn{D^{(1)}} and diffusion
#' \eqn{D^{(2)}} directly from measured data.
#'
#' \deqn{ \dot{X}(t) = D^{(1)}(X(t),t) + \sqrt{D^{(2)}(X(t),t)}\,\Gamma(t)\quad
#' \mathrm{with} } \deqn{ D^{(n)}(x,t) = \lim_{\tau \rightarrow 0}
#' \frac{1}{\tau} M^{(n)}(x,t,\tau)\quad \mathrm{and} } \deqn{ M^{(n)}(x,t,\tau)
#' = \frac{1}{n!} \langle (X(t+\tau) - x)^n \rangle |_{X(t) = x} }
#'
#' The Langevin equation should be interpreted in the way that for every time
#' \eqn{t_i} where the system meets an arbitrary but fixed point \eqn{x} in
#' phase space, \eqn{X(t_i+\tau)} is defined by the deterministic function
#' \eqn{D^{(1)}(x)} and the stochastic function
#' \eqn{\sqrt{D^{(2)}(x)}\Gamma(t_i)}. Both, \eqn{D^{(1)}(x)} and
#' \eqn{D^{(2)}(x)} are constant for fixed \eqn{x}.
#'
#' One can integrate drift and diffusion numerically over small intervals. If
#' the system is at time \eqn{t} in the state \eqn{x = X(t)} the drift can be
#' calculated for small \eqn{\tau} by averaging over the difference of the
#' system state at \eqn{t+\tau} and the state at \eqn{t}. The average has to be
#' taken over the whole ensemble or in the stationary case over all \eqn{t =
#' t_i} with \eqn{X(t_i) = x}. Diffusion can be calculated analogously.
#'
#'
#' @name Langevin-package
#' @aliases Langevin-package
#' @docType package
#' @author Philip Rinn
#' @references
#'
#' \bold{A review of the Langevin method can be found at:}
#'
#' Friedrich, R.; et al. (2011) \emph{Approaching Complexity by Stochastic
#' Methods: From Biological Systems to Turbulence}. Physics Reports, 506(5), 87–162.
#'
#' \bold{For further reading:}
#'
#' Risken, H. (1996) \emph{The Fokker-Planck equation}. Springer.
#'
#' Siegert, S.; et al. (1998) \emph{Analysis of data sets of stochastic
#' systems}. Phys. Lett. A.
#'
#' Friedrich, R.; et al. (2000) \emph{Extracting model equations from
#' experimental data}. Phys. Lett. A.
NULL
