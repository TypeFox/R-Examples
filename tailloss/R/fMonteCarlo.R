#' Monte Carlo Simulations.
#'
#' Function to estimate the total losses via the Monte Carlo simulations.
#'
#' @param ELT Data frame containing two numeric columns. The column \code{Loss} contains the expected losses from each single occurrence of event. The column \code{Rate} contains the arrival rates of a single occurrence of event.
#' @param s Scalar or numeric vector containing the total losses of interest. 
#' @param t Scalar representing the time period of interest. The default value is \code{t} = 1.
#' @param theta Scalar containing information about the variance of the Gamma distribution: \eqn{sd[X] = x * }\code{theta}. The default value is \code{theta} = 0: the loss associated to an event is considered as a constant.
#' @param cap Scalar representing the financial cap on losses for a single event, i.e. the maximum possible loss caused by a single event. The default value is \code{cap} = Inf.
#' @param nsim Integer representing the number of Monte Carlo simulations. The default value is \code{nsim} = 10e3.
#' @param verbose Logical, if \code{TRUE} returns 95\% CB and raw sample. The default is \code{verbose = FALSE}.

#' @return If \code{verbose = FALSE} the function returns a numeric matrix, containing in the first column the pre-specified losses \code{s}, and the estimated exceedance probabilities in the second column.
#' If \code{verbose = TRUE} the function returns a numeric matrix  containing four columns. The first column contains the losses \code{s}, the second column contains the estimated exceedance probabilities, the other columns contain the 95\% confidence bands. The attributes of this matrix are a vector \code{simS} containing the simulated losses.

#' @keywords montecarlo
#' @export
#' @examples
#' data(UShurricane)
#'
#' # Compress the table to millions of dollars
#'
#' USh.m <- compressELT(ELT(UShurricane), digits = -6)
#' EPC.MonteCarlo <- fMonteCarlo(USh.m, s = 1:40, verbose = TRUE)
#' EPC.MonteCarlo
#' par(mfrow = c(1, 2))
#' plot(EPC.MonteCarlo[, 1:2], type = "l", ylim = c(0, 1))
#' matlines(EPC.MonteCarlo[, -2], ylim = c(0, 1), lty = 2, col = 1)
#' # Assuming the losses follow a Gamma with E[X] = x, and Var[X] = 2 * x and cap = 5m
#' EPC.MonteCarlo.Gamma <- fMonteCarlo(USh.m, s = 1:40, theta = 2, cap = 5, verbose = TRUE)
#' EPC.MonteCarlo.Gamma
#' plot(EPC.MonteCarlo.Gamma[, 1:2], type = "l", ylim = c(0, 1))
#' matlines(EPC.MonteCarlo.Gamma[, -2], ylim = c(0,1), lty = 2, col = 1)
#' # Compare the two results:
#' par(mfrow = c(1, 1))
#' plot(EPC.MonteCarlo[, 1:2], type = "l", main = "Exceedance Probability Curve", 
#' ylim = c(0, 1))
#' lines(EPC.MonteCarlo.Gamma[, 1:2], col = 2, lty = 2)
#' legend("topright", c("Dirac Delta", expression(paste("Gamma(", 
#' alpha[i] == 1 / theta^2, ", ", beta[i] ==1 / (x[i] * theta^2), ")", " cap =", 5))),  
#' lwd = 2, lty = 1:2, col = 1:2)

fMonteCarlo <- function(ELT, s, t = 1, theta = 0, cap = Inf, nsim = 10e3, verbose = FALSE) {

		stopifnot(inherits(ELT, "ELT"),
            s >= 0, 
            is.scalar(t), t > 0,
            is.scalar(theta), theta >= 0, 
            is.scalar(cap), cap > 0,
            is.scalar(nsim), nsim > 0, nsim == floor(nsim),
            is.logical(verbose))

  ## generate empirical distribution of losses

  m <- length(ELT$Loss)
  lambda <- sum(ELT$Rate)
  N <- rpois(nsim, lambda * t)

  ##  Branch on theta for speed

  if (theta == 0) {
  	  	
  	 if (!is.infinite(cap))
  	 	ELT$Loss <- pmin(cap, ELT$Loss)
  	
    losses <- sapply(1L:nsim, function(j) {
      ii <- sample.int(m, N[j], replace = TRUE, prob = ELT$Rate)
      xx <- ELT$Loss[ii]
      sum(xx)
    })

  } else {

    alpha <- 1 / theta^2
    beta <- alpha / ELT$Loss
    
    losses <- sapply(1L:nsim, function(j) {
    		ii <- sample.int(m, N[j], replace = TRUE, prob = ELT$Rate)
    		xx <- rgamma(N[j], shape = alpha, rate = beta[ii])
    		if(!is.infinite(cap)) 
    				xx <- pmin(cap, xx)
    		sum(xx)
    		})
      
  }

  ## compute the exceedance probabilities and return

  pp <- colMeans(outer(losses, s, ">="))
  robj <- cbind(s = s, "Pr[S>=s]" = pmin(1, pp))

  if (verbose) {

    ## attach 95% confidence bounds, pp is prob of exceedance

    alpha <- 0.05
    epsilon <- sqrt(log(2 / alpha) / (2 * nsim))
    L <- pmax(0, (1 - pp) - epsilon)
    U <- pmin(1, (1 - pp) + epsilon)
    robj <- cbind(robj, L95 = 1 - U, U95 = 1 - L)
    attr(robj, "rsam") <- losses
  }

  robj
}