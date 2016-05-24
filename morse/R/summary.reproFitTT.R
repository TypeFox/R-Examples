#' Summary for reproFitTT objects
#'
#' The summary shows the quantiles of priors and posteriors on parameters
#' and the quantiles of the posterior on the ECx.
#'
#' @param object an object of class \code{reproFitTT}
#' @param quiet when \code{FALSE}, prints summary on standard output
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @return The function returns a list with the following fields:
#' \item{Qpriors}{quantiles for the model's prior}
#' \item{Qposteriors}{quantiles for the model's posteriors}
#' \item{QECx}{quantiles for ECx values}
#'
#' @seealso reproFitTT
#'
#' @examples
#' # (1) Load the data
#' data(cadmium1)
#'
#' # (2) Create a reproData object
#' cadmium1 <- reproData(cadmium1)
#'
#' \dontrun{
#' # (3) Run the reproFitTT function with the log-logistic
#' # model
#' out <- reproFitTT(dat, ecx = c(5, 10, 15, 20, 30, 50, 80),
#' quiet = TRUE)
#'
#' # (4) summarize the reproFitTT object
#' summary(out)
#' }
#'
#' @keywords summary
#' 
#' @importFrom stats qnorm qunif
#' 
#' @export
summary.reproFitTT <- function(object, quiet = FALSE, ...) {

  # quantiles of priors parameters
  n.iter <- object$n.iter$end - object$n.iter$start

  # b
  log10b <- qunif(p = c(0.5, 0.025, 0.975),
                  min = object$jags.data$log10bmin,
                  max = object$jags.data$log10bmax)

  b <- 10^log10b

  # d
  d <- qnorm(p = c(0.5, 0.025, 0.975),
             mean = object$jags.data$meand,
             sd = 1 / sqrt(object$jags.data$taud))

  # e
  log10e <- qnorm(p = c(0.5, 0.025, 0.975),
                  mean = object$jags.data$meanlog10e,
                  sd = 1 / sqrt(object$jags.data$taulog10e))

  e <- 10^log10e

  if (object$model.label == "P") {
    res <- rbind(b, d, e)
  }
  if (object$model.label == "GP") {
    # omega
    log10omega <- qunif(p = c(0.5, 0.025, 0.975),
                        min = object$jags.data$log10omegamin,
                        max = object$jags.data$log10omegamax)

    omega <- 10^log10omega

    res <- rbind(b, d, e, omega)
  }

  ans1 <-  round(data.frame(res), digits = 3)
  colnames(ans1) <- c("50%", "2.5%", "97.5%")

  # quantiles of estimated model parameters
  ans2 <- round(object$estim.par, digits = 3)
  colnames(ans2) <- c("50%", "2.5%", "97.5%")

  # estimated ECx and their CIs 95%
  ans3 <- round(object$estim.ECx, digits = 3)
  colnames(ans3) <- c("50%", "2.5%", "97.5%")

  if (! quiet) {
    cat("Summary: \n\n")
    if (object$model.label == "GP")
      cat("The ", object$det.part, " model with a Gamma Poisson stochastic part was used !\n\n")
    if(object$model.label == "P")
      cat("The ", object$det.part, " model with a Poisson stochastic part was used !\n\n")
    cat("Priors on parameters (quantiles):\n\n")
    print(ans1)
    cat("\nPosterior of the parameters (quantiles):\n\n")
    print(ans2)
    cat("\nPosterior of the ECx (quantiles):\n\n")
    print(ans3)
  }

  invisible(list(Qpriors = ans1,
                 Qposteriors = ans2,
                 QECx = ans3))
}

