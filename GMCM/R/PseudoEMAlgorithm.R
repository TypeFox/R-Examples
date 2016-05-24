#' EM-like algorithm for the GMCM
#'
#' An fast and modified implementation of the Li et. al. (2011) EM-like
#' algorithm for estimating the maximizing parameters of the GMCM-likelihood
#' function.
#'
#' When either \code{"absGMCM"} or \code{"absLi"} are used, the parameters
#' corresponding to the biggest observed likelihood is returned. This is not
#' necessarily the last iteration.
#'
#' @param x A matrix of observations where rows corresponds to features and
#' columns to experiments.
#' @param theta A list of parameters formatted as described in
#' \code{\link{rtheta}}.
#' @param eps The maximum difference required to achieve convergence.
#' @param max.ite The maximum number of iterations.
#' @param verbose Logical. Set to \code{TRUE} to increase verbosity.
#' @param trace.theta Logical. If \code{TRUE}, a trace of the estimated thetas
#' are returned.
#' @param meta.special.case Logical. If \code{TRUE}, the estimators used are for
#' the special case GMCM of Li et. al. (2011).
#' @param convergence.criterion Character. Sets the convergence criterion.  If
#' \code{"absGMCM"} the absolute value of difference in GMCM is used. If
#' \code{"GMCM"} the difference in GMCM-likelihoods are used as convergence
#' criterion. If \code{"GMM"}, the guaranteed non-decreasing difference of
#' GMM-likelihoods are used. If \code{"Li"}, the convergence criterion used by
#' Li et. al. (2011) is used. If \code{"absLi"}, the absolute values of the Li
#' et. al. criterion.
#' @return A list of 3 or 4 is returned depending on the value of
#' \code{trace.theta} \item{theta}{A list containing the final parameter
#' estimate in the format of \code{\link{rtheta}}} \item{loglik.tr}{A matrix
#' with different log-likelihood traces in each row.} \item{kappa}{A matrix
#' where the (i,j)'th entry is the probability that \code{x[i, ]} belongs to
#' the j'th component. Usually the returned value of \code{EStep}.}
#' \item{theta.tr}{A list of each obtained parameter estimates in the format of
#' \code{\link{rtheta}}}
#' @note The algorithm is highly sensitive to the starting parameters which
#' therefore should be carefully chosen.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @seealso \code{\link{rtheta}}, \code{\link{EMAlgorithm}}
#' @references
#'   Li, Q., Brown, J. B. J. B., Huang, H., & Bickel, P. J. (2011).
#'   Measuring reproducibility of high-throughput experiments. The Annals of
#'   Applied Statistics, 5(3), 1752-1779. doi:10.1214/11-AOAS466
#' @examples
#' set.seed(1)
#'
#' # Choosing the true parameters and simulating data
#' true.par <- c(0.8, 3, 1.5, 0.4)
#' data <- SimulateGMCMData(n = 1000, par = true.par, d = 2)
#' uhat <- Uhat(data$u)  # Observed ranks
#'
#' # Plot of latent and observed data colour coded by the true component
#' par(mfrow = c(1,2))
#' plot(data$z, main = "Latent data", cex = 0.6,
#'      xlab = "z (Experiment 1)", ylab = "z (Experiment 2)",
#'      col = c("red","blue")[data$K])
#' plot(uhat, main = "Observed data", cex = 0.6,
#'      xlab = "u (Experiment 1)", ylab = "u (Experiment 2)",
#'      col = c("red","blue")[data$K])
#'
#' # Fit the model using the Pseudo EM algorithm
#' init.par <- c(0.5, 1, 1, 0.5)
#' res <- GMCM:::PseudoEMAlgorithm(uhat, meta2full(init.par, d = 2),
#'                                 verbose = TRUE,
#'                                 convergence.criterion = "absGMCM",
#'                                 eps = 1e-4,
#'                                 trace.theta = FALSE,
#'                                 meta.special.case = TRUE)
#'
#' # Compute posterior cluster probabilities
#' IDRs <- get.IDR(uhat, par = full2meta(res$theta))
#'
#' # Plot of observed data colour coded by the MAP estimate
#' plot(res$loglik[3,], main = "Loglikelihood trace", type = "l",
#'      ylab = "log GMCM likelihood")
#' abline(v = which.max(res$loglik[3,])) # Chosen MLE
#' plot(uhat, main = "Clustering\nIDR < 0.05", xlab = "", ylab = "", cex = 0.6,
#'      col = c("Red","Blue")[IDRs$Khat])
#'
#' # View parameters
#' rbind(init.par, true.par, estimate = full2meta(res$theta))
#'
#' # Confusion matrix
#' table("Khat" = IDRs$Khat, "K" = data$K)
#' @keywords internal
PseudoEMAlgorithm <- function (x, theta,
                               eps = 1e-4,
                               max.ite = 1e3,
                               verbose = FALSE,
                               trace.theta = FALSE,
                               meta.special.case = FALSE,
                               convergence.criterion =
                                 c("absGMCM", "GMCM", "GMM", "Li", "absLi")) {
  cc <- match.arg(convergence.criterion)

  u <- Uhat(x)
  z <- qgmm.marginal(u, theta)
  loglik.tr           <- matrix(NA, 3, max.ite)
  loglik.tr[2, 1]     <- dgmm.loglik(theta, z)
  loglik.tr[3, 1]     <- dgmcm.loglik(theta, u)
  rownames(loglik.tr) <- c("gmm.pre", "gmm.post", "gmcm")
  theta.tr            <- list(theta)

  for (k in 2:max.ite) {
    z <- qgmm.marginal(u, theta)
    loglik.pre <- dgmm.loglik(theta, z)    # Compute loglik pre EM step
    kappa <- EStep(x = z, theta = theta)
    if (any(colSums(kappa) == 0)) {
      stop("No observations are estimated to be from component(s): ",
           paste(which(colSums(kappa) == 0), collapse = " and "), ". ",
           "All posterior probabilities are zero. ",
           "Try another start estimate or fewer components.")
    }
    theta <- MStep(x = z, kappa = kappa,
                   meta.special.case = meta.special.case)
    loglik.post      <- dgmm.loglik(theta, z)   # Compute loglik post EM step
    loglik.tr[1:2,k] <- c(loglik.pre, loglik.post)  # Tracking of pre/post theta
    loglik.tr[3,  k] <- dgmcm.loglik(theta, u)
    theta.tr[[k]] <- theta

    delta <-
      switch(cc,
             "GMCM"    = loglik.tr[3, k] - loglik.tr[3, k-1],
             "absGMCM" = abs(loglik.tr[3, k] - loglik.tr[3, k-1]),
             "GMM"     = loglik.post - loglik.pre,
             "Li"      = loglik.tr[2, k]- loglik.tr[3, k-1],
             "absLi"   = abs(loglik.tr[2, k]- loglik.tr[3, k-1]))

    if (verbose) {
      cat(k,
          "| delta =", sprintf("%.6f", delta),
          "| gmm =", sprintf("%.2f", round(loglik.post,2)),
          "| gmcm = ", sprintf("%.2f", round(loglik.tr[3,k],2)),"\n")
      flush.console()
    }

    if (delta < eps) {
      break
    }
  }

  if (k == max.ite)  # Warn if maximun interations reached.
    warning(paste("Max iterations (", max.ite, ") reached.", sep = ""))

  if (cc %in% c("absGMCM", "absLi")) {
    k <- which.max(loglik.tr[3, ])
    theta <- theta.tr[[k]]
  }
  loglik.tr <- loglik.tr[, !apply(is.na(loglik.tr), 2, all)]
  res <- list(theta = theta, loglik.tr = loglik.tr,
              kappa = kappa, theta.tr  = theta.tr)
  if (!trace.theta)
    res <- res[-4]
  return(res)
}
