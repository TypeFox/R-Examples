#' Reproducibility/meta analysis using GMCMs
#'
#' This function performs reproducibility (or meta) analysis using GMCMs.
#' It features various optimization routines to identify the maximum likelihood
#' estimate of the special Gaussian mixture copula model proposed by
#' Li et. al. (2011).
#'
#' The \code{"L-BFGS-B"} method does not perform a transformation of the
#' parameters.
#'
#' @aliases fit.meta.gmcm
#' @param u An \code{n} by \code{d} matrix of test statistics. Rows correspond
#'   to features and columns to experiments. Large values are assumed to be
#'   indicative of reproducible genes.
#' @param init.par A 4-dimensional vector of the initial parameters where,
#'   \code{init.par[1]} is the mixture proportion of spurious signals,
#'   \code{init.par[2]} is the mean, \code{init.par[3]} is the standard
#'   deviation, \code{init.par[4]} is the correlation.
#' @param method A character vector of length \eqn{1}{1}. The optimization
#'   method used. Should be either \code{"NM"}, \code{"SANN"}, \code{"L-BFGS"},
#'   \code{"L-BFGS-B"}, or \code{"PEM"} which are abbreviations of Nelder-Mead,
#'   Simulated Annealing, limited-memory quasi-Newton method, limited-memory
#'   quasi-Newton method with box constraints, and the pseudo EM algorithm,
#'   respectively. Default is \code{"NM"}. See \code{\link{optim}} for further
#'   details.
#' @param max.ite The maximum number of iterations.  If the \code{method} is
#'   \code{"SANN"} this is the number of iterations as there is no other
#'   stopping criterion. (See \code{\link{optim}})
#' @param verbose Logical. If \code{TRUE}, the log-likelihood values are
#'   printed.
#' @param positive.rho \code{logical}. If \code{TRUE}, the correlation parameter
#'   is restricted to be positive.
#' @param trace.theta \code{logical}. Extra convergence information is appended
#'   as a list to the output returned if \code{TRUE}. The exact behavior is
#'   dependent on the value of \code{method}. If \code{method} equals
#'   \code{"PEM"}, the argument is passed to \code{trace.theta} in
#'   \code{\link{PseudoEMAlgorithm}}. Otherwise it is passed to the control
#'   argument \code{trace} in \code{\link{optim}}.
#' @param \dots Arguments passed to the \code{control}-list in
#'   \code{\link{optim}} or \code{\link{PseudoEMAlgorithm}} if \code{method} is
#'   \code{"PEM"}.
#' @return A vector \code{par} of length 4 of the fitted parameters where
#'   \code{par[1]} is the probability of being from the first (or null)
#'   component, \code{par[2]} is the mean, \code{par[3]} is the standard
#'   deviation, and \code{par[4]} is the correlation.
#' @note Simulated annealing is strongly dependent on the initial values and
#'   the cooling scheme.
#'
#'   See \code{\link{optim}} for further details.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @seealso \code{\link{optim}}
#' @references
#'   Li, Q., Brown, J. B. J. B., Huang, H., & Bickel, P. J. (2011).
#'   Measuring reproducibility of high-throughput experiments. The Annals of
#'   Applied Statistics, 5(3), 1752-1779. doi:10.1214/11-AOAS466
#' @examples
#' set.seed(1)
#'
#' # True parameters
#' true.par <- c(0.9, 2, 0.7, 0.6)
#' # Simulation of data from the GMCM model
#' data <- SimulateGMCMData(n = 1000, par = true.par)
#' uhat <- Uhat(data$u) # Ranked observed data
#'
#' init.par <- c(0.5, 1, 0.5, 0.9)  # Initial parameters
#'
#' # Optimization with Nelder-Mead
#' nm.par   <- fit.meta.GMCM(uhat, init.par = init.par, method = "NM")
#'
#' \dontrun{
#' # Comparison with other optimization methods
#' # Optimization with simulated annealing
#' sann.par <- fit.meta.GMCM(uhat, init.par = init.par, method = "SANN",
#'                           max.ite = 3000, temp = 1)
#' # Optimization with the Pseudo EM algorithm
#' pem.par  <- fit.meta.GMCM(uhat, init.par = init.par, method = "PEM")
#'
#' # The estimates agree nicely
#' rbind("True" = true.par, "Start" = init.par,
#'       "NM" = nm.par, "SANN" = sann.par, "PEM" = pem.par)
#' }
#'
#' # Get estimated cluster
#' Khat <- get.IDR(x = uhat, par = nm.par)$Khat
#' plot(uhat, col = Khat, main = "Clustering\nIDR < 0.05")
#' @export
fit.meta.GMCM <- function(u,
                          init.par,
                          method = c("NM", "SANN", "L-BFGS", "L-BFGS-B", "PEM"),
                          max.ite = 1000,
                          verbose = TRUE,
                          positive.rho = TRUE,
                          trace.theta = FALSE,
                          ...)
{
  d <- ncol(u)

  # Note, Uhat is idempotent. Hence, already ranked data will not change
  u <- Uhat(u)

  switch(match.arg(method),
         "NM" = {  # Fitted using Nelder-Mead (The amoeba method)
           fit <- optim(inv.tt(init.par, d = d, positive.rho = positive.rho),
                        meta.gmcm.loglik, u = u,
                        positive.rho = positive.rho,
                        rescale = TRUE, method = "Nelder-Mead",
                        control = list(maxit = max.ite,
                                       trace = verbose,
                                       fnscale = -1, ...))
           fitted.par <- tt(fit$par, d = d, positive.rho = positive.rho)
         },
         "SANN" = {  # Fitting using Simulated Annealing
           fit <- optim(inv.tt(init.par, d = d, positive.rho = positive.rho),
                        meta.gmcm.loglik,  u = u, positive.rho = positive.rho,
                        rescale = TRUE, method = "SANN",
                        control = list(maxit = max.ite,
                                       trace = verbose,
                                       fnscale = -1, ...))
           fitted.par <- tt(fit$par, d = d, positive.rho = positive.rho)
         },
         "L-BFGS" = { # Fit via L-BFGS-B (Limited-memory quasi-Newton method)
           fit <- optim(inv.tt(init.par, d = d,
                               positive.rho = positive.rho),
                        meta.gmcm.loglik, u = u, positive.rho = positive.rho,
                        rescale = TRUE, method = "L-BFGS-B",
                        control = list(maxit = max.ite,
                                       trace = verbose,
                                       fnscale = -1, ...))
           fitted.par <- tt(fit$par, d = d, positive.rho = positive.rho)
         },
         "L-BFGS-B" = {  # Fitting using L-BFGS-B !! NOT RESCALED !!
           fit <- optim(init.par,
                        meta.gmcm.loglik, u = u,
                        positive.rho = positive.rho,
                        rescale = FALSE, method = "L-BFGS-B",
                        upper = c(1,Inf,Inf,1),
                        lower = c(0,0,0, ifelse(positive.rho, 0, -1/(d-1))),
                        control = list(maxit = max.ite,
                                       trace = verbose,
                                       fnscale = -1, ...))
           fitted.par <- fit$par
         },
         "PEM" = {  # Fitting using the Li Pseudo EM Algorithm
           fit <- PseudoEMAlgorithm(u, theta = meta2full(init.par, d = d),
                                    max.ite = max.ite,
                                    meta.special.case = TRUE,
                                    trace.theta = trace.theta,
                                    verbose = verbose,
                                    ...)
           fitted.par <- full2meta(fit$theta)
         })
  names(fitted.par) <- c("pie1", "mu", "sigma", "rho")
  if (trace.theta)
    fitted.par <- list(fitted.par, fit)
  return(fitted.par)
}
