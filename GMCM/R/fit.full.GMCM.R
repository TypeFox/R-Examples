#' Unsupervised clustering using a general GMCM
#'
#' Perform unsupervised clustering using various optimization procedures to find
#' the maximum likelihood estimate of the general Gaussian mixture copula
#' model by Tewari et al. (2011).
#'
#' The \code{"L-BFGS-B"} method does not perform a transformation of the
#' parameters and uses box constraints as implemented in \code{optim}. \cr
#' Note that the many parameter configurations are poorly estimable or directly
#' unidentifiable.
#'
#' @aliases fit.full.gmcm
#' @param u An \code{n} by \code{d} matrix of ranked and scaled test statistics.
#'   Rows
#'   correspond to observations and columns to the dimensions of the variables.
#' @param m The number of components to be fitted.
#' @param theta A list of parameters as defined in \code{\link{rtheta}}. If
#'   \code{theta} is not provided, then heuristic starting values are chosen
#'   using the k-means algorithm.
#' @param method A character vector of length \eqn{1}{1}. The optimization
#'   method used. Should be either \code{"NM"}, \code{"SANN"}, \code{"L-BFGS"},
#'   \code{"L-BFGS-B"}, or \code{"PEM"} which are the Nelder-Mead, Simulated
#'   Annealing, limited-memory quasi-Newton method, limited-memory quasi-Newton
#'   method with box constraints, and the pseudo EM algorithm, respectively.
#'   Default is \code{"NM"}. See \code{\link{optim}} for further details.
#' @param max.ite The maximum number of iterations. If the \code{method} is
#'   \code{"SANN"} this is the number of iterations as there is no other
#'   stopping criterion. (See \code{\link{optim}})
#' @param verbose Logical. If \code{TRUE}, a trace of the parameter estimates
#'   is made.
#' @param \dots Arguments passed to the \code{control}-list in
#'   \code{\link{optim}} or \code{\link{PseudoEMAlgorithm}} if the \code{method}
#'   is \code{"PEM"}.
#' @return A list of parameters formatted as described in \code{\link{rtheta}}.
#' @note All the optimization procedures are strongly dependent on the initial
#'   values and the cooling scheme. Therefore it is advisable to apply multiple
#'   different initial parameters and select the best fit.
#'
#'   The \code{\link{choose.theta}} itself chooses random a initialization.
#'   Hence, the output when \code{theta} is not directly supplied can vary.
#'
#'   See \code{\link{optim}} for further details.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @seealso \code{\link{optim}}, \code{\link{get.prob}}
#' @references
#'   Li, Q., Brown, J. B. J. B., Huang, H., & Bickel, P. J. (2011).
#'   Measuring reproducibility of high-throughput experiments. The Annals of
#'   Applied Statistics, 5(3), 1752-1779. doi:10.1214/11-AOAS466
#'
#'   Tewari, A., Giering, M. J., & Raghunathan, A. (2011). Parametric
#'   Characterization of Multimodal Distributions with Non-gaussian Modes. 2011
#'   IEEE 11th International Conference on Data Mining Workshops, 286-292.
#' doi:10.1109/ICDMW.2011.135
#' @examples
#' set.seed(17)
#' sim <- SimulateGMCMData(n = 1000, m = 3, d = 2)
#'
#' # Plotting simulated data
#' par(mfrow = c(1,2))
#' plot(sim$z, col = rainbow(3)[sim$K], main = "Latent process")
#' plot(sim$u, col = rainbow(3)[sim$K], main = "GMCM process")
#'
#' # Observed data
#' uhat <- Uhat(sim$u)
#'
#' # The model should be fitted multiple times using different starting estimates
#' start.theta <- choose.theta(uhat, m = 3)  # Random starting estimate
#' res <- fit.full.GMCM(u = uhat, theta = start.theta,
#'                      method = "NM", max.ite = 3000,
#'                      reltol = 1e-2, trace = TRUE)  # Note, 1e-2 is too big
#'
#' # Confusion matrix
#' Khat <- apply(get.prob(uhat, theta = res), 1, which.max)
#' table("Khat" = Khat, "K" = sim$K)  # Note, some components have been swapped
#'
#' # Simulation from GMCM with the fitted parameters
#' simfit <- SimulateGMCMData(n = 1000, theta = res)
#'
#' # As seen, the underlying latent process is hard to estimate.
#' # The clustering, however, is very good.
#' par(mfrow = c(2,2))
#' plot(simfit$z, col = simfit$K, main = "Model check 1\nSimulated GMM")
#' plot(simfit$u, col = simfit$K, main = "Model check 2\nSimulated GMCM")
#' plot(sim$u, col = Khat, main = "MAP clustering")
#' @export
fit.full.GMCM <- function (u,
                           m,
                           theta = choose.theta(u, m),
                           method = c("NM", "SANN", "L-BFGS", "L-BFGS-B", "PEM"),
                           max.ite = 1000,
                           verbose = TRUE,
                           ...) {
  # Note, Uhat is idempotent. Hence, already ranked data will not change
  u <- Uhat(u)

  method <- gsub("NM", "Nelder-Mead", match.arg(method))

  if (missing(m) & missing(theta)) {
    stop("m is not supplied.")
  }

  if (method != "PEM") {

    gmcm.loglik <- function (par, u, m) { # Defining objective function
      #cat("par=", par, "\n")  # FOR DEBUGGING
      theta <- vector2theta(par, d = ncol(u), m = m)
      theta$pie <- theta$pie/sum(theta$pie)
      #cat("theta=", unlist(theta$sigma), "\n")  # FOR DEBUGGING
      loglik <- dgmcm.loglik(theta = theta, u)

      return(loglik)
    }

    par <- theta2vector(theta)
    fit <- optim(par, gmcm.loglik, u = u, m = theta$m,
                 control = list(maxit = max.ite,
                                fnscale = -1, trace = verbose, ...),
                 method = method)
    theta <- vector2theta(fit$par, d = theta$d, m = theta$m)
    theta$pie <- theta$pie/sum(theta$pie)

    return(theta)

  } else {

    fit <- PseudoEMAlgorithm(x = u,
                             theta = theta,
                             max.ite = max.ite,
                             verbose = verbose,
                             meta.special.case = FALSE,
                             ...)

    return(fit$theta)
  }
}
