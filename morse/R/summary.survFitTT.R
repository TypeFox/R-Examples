#' Summary for survFitTT objects
#'
#' The summary shows the quantiles of priors and posteriors on parameters
#' and the quantiles of the posterior on the LCx.
#'
#' @param object an object of class \code{survFitTT}
#' @param quiet when \code{FALSE}, prints summary on standard output
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @return The function returns a list with the following fields:
#' \item{Qpriors}{quantiles for the model's prior}
#' \item{Qposteriors}{quantiles for the model's posteriors}
#' \item{QLCx}{quantiles for LCx values}
#'
#' @seealso survFitTT
#'
#' @examples
#' # (1) Load the data
#' data(cadmium1)
#'
#' # (2) Create a survData object
#' cadmium1 <- survData(cadmium1)
#'
#' \dontrun{
#' # (3) Run the survFitTT function with the log-logistic
#' # binomial model
#' out <- survFitTT(cadmium1, lcx = c(5, 10, 15, 20, 30, 50, 80),
#'                  quiet = TRUE)
#'
#' # (4) summarize the survFitTT object
#' summary(out)
#' }
#'
#' @keywords summary
#'
#' @importFrom stats qnorm qunif
#' 
#' @export
summary.survFitTT <- function(object, quiet = FALSE, ...) {

  # quantiles of priors parameters
  n.iter <- object$n.iter$end - object$n.iter$start

  # b
  log10b <- qunif(p = c(0.5, 0.025, 0.975),
                  min = object$jags.data$log10bmin,
                  max = object$jags.data$log10bmax)

  b <- 10^log10b

  # e
  log10e <- qnorm(p = c(0.5, 0.025, 0.975),
                  mean = object$jags.data$meanlog10e,
                  sd = 1 / sqrt(object$jags.data$taulog10e))

  e <- 10^log10e

  # d
  if (object$det.part == "loglogisticbinom_3") {

    d <- qunif(p = c(0.5, 0.025, 0.975),
               min = object$jags.data$dmin,
               max = object$jags.data$dmax)

    res <- rbind(b, d, e)
  } else {
    res <- rbind(b, e)
  }

  ans1 <- round(data.frame(res), digits = 3)
  colnames(ans1) <- c("50%", "2.5%", "97.5%")

  # quantiles of estimated model parameters
  ans2 <- round(object$estim.par, digits = 3)
  colnames(ans2) <- c("50%", "2.5%", "97.5%")

  # estimated ECx and their CIs 95%
  ans3 <- round(object$estim.LCx, digits = 3)
  colnames(ans3) <- c("50%", "2.5%", "97.5%")

  # print
  if (! quiet) {
    cat("Summary: \n\n")
    cat("The ", object$det.part, " model with a binomial stochastic part was used !\n\n")
    cat("Priors on parameters (quantiles):\n\n")
    print(ans1)
    cat("\nPosterior of the parameters (quantiles):\n\n")
    print(ans2)
    cat("\nPosterior of the LCx (quantiles):\n\n")
    print(ans3)
  }

  invisible(list(Qpriors = ans1,
                 Qpost = ans2,
                 QLCx = ans3))
}

