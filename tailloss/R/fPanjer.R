#' Panjer Recursion.
#'
#' Function to calculate the total losses via the Panjer recursion.
#'
#' @param ELT Data frame containing two numeric columns. The column \code{Loss} contains the expected losses from each single occurrence of event. The column \code{Rate} contains the arrival rates of a single occurrence of event.
#' @param s Scalar or numeric vector containing the total losses of interest. 
#' @param t Scalar representing the time period of interest. The default value is \code{t} = 1.
#' @param theta Scalar containing information about the variance of the Gamma distribution: \eqn{sd[X] = x * }\code{theta}. The default value is \code{theta} = 0: the loss associated to an event is considered as a constant.
#' @param cap Scalar representing the financial cap on losses for a single event, i.e. the maximum possible loss caused by a single event. The default value is \code{cap} = Inf.
#' @param nq Scalar, number of quantiles added when \code{theta > 0}
#' @param verbose A logical, if \code{TRUE} gives the entire distribution up to the maximum value of \code{s}. If \code{FALSE} gives only the results for the specified values of \code{s}. The default is \code{verbose = FALSE}.

#'  @return A numeric matrix containing the pre-specified losses \code{s}  in the first column and the exceedance probabilities in the second column.

#' @references Panjer, H.H. (1980), `The aggregate claims distribution and stop-loss reinsurance,' \emph{Transactions of the Society of Actuaries}, 32, 523-545.

#' @keywords Panjer
#' @export
#' @examples
#' data(UShurricane)
#'
#' # Compress the table to millions of dollars
#'
#' USh.m <- compressELT(ELT(UShurricane), digits = -6)
#'
#' EPC.Panjer <- fPanjer(USh.m, s = 1:40, verbose = TRUE)
#' EPC.Panjer
#' plot(EPC.Panjer, type = "l", ylim = c(0,1))
#' # Assuming the losses follow a Gamma with E[X] = x, and Var[X] = 2 * x and cap = 5m
#'
#' EPC.Panjer.Gamma <- fPanjer(USh.m, s = 1:40, theta = 2, cap = 5, verbose = TRUE)
#' EPC.Panjer.Gamma
#' plot(EPC.Panjer.Gamma, type = "l", ylim = c(0,1))
#'
#' # Compare the two results:
#'
#' plot(EPC.Panjer, type = "l", main = 'Exceedance Probability Curve', 
#' ylim = c(0, 1))
#' lines(EPC.Panjer.Gamma, col = 2, lty = 2)
#' legend("topright", c("Dirac Delta", expression(paste("Gamma(", 
#' alpha[i] == 1 / theta^2, ", ", beta[i] ==1 / (x[i] * theta^2), ")", " cap =", 5))),  
#' lwd = 2, lty = 1:2, col = 1:2)

fPanjer <- function(ELT, s, t = 1, theta = 0, cap = Inf, nq = 10, verbose = FALSE){ # Panjer recursion

	stopifnot(inherits(ELT, "ELT"),
            s >= 0, 
            is.scalar(t), t > 0,
            is.scalar(theta), theta >= 0,
            is.scalar(cap), cap > 0,
            is.logical(verbose))

  if (theta == 0) {
  	  	
  	 if (!is.infinite(cap)){
  	 	ELT$Loss <- pmin(cap, ELT$Loss)
  	 	ELT <- compressELT(ELT, digits = 0)
  	 	}

  } else {   ## handle case of theta > 0 as new ELT

    stopifnot(is.scalar(nq), nq > 0, nq == floor(nq))
    
    	alpha <- 1 / theta^2
    	beta <- alpha / ELT$Loss
    	qq <- 1 / nq
    	qq <- seq(from = qq / 2, by = qq, length = nq)
    	
    	loss <- qgamma(qq, shape = alpha, rate = rep(beta, each = nq))

    	if (!is.infinite(cap))
    		loss <- pmin(loss, cap)
      		
		rate <- rep(ELT$Rate / nq, each = nq)
    	ELT <- compressELT(ELT(Rate = rate, Loss = loss), digits = 0)
   
  }
  
  ## Panjer recursion

  m <- max(ELT$Loss)
  lambda <- sum(ELT$Rate)
  alpha <- rep(0, m)
  alpha[ELT$Loss] <- ELT$Rate / lambda # alpha[j] = Pr(X = j), j = 1, 2, ..., m
  P <- rep(NA, max(s) + 1) # P[n+1] = P_n, n = 0, 1, ...

  lambda <- lambda * t
  P[1] <- exp(-lambda) # n = 0
  for (n in 1L:max(s)) {
      j <- 1L:min(m, n)
      P[n+1L] <- lambda * sum(j * alpha[j] * P[n-j+1L]) / n
  }

  ## extract the values for s

  pp <- c(1, 1 - cumsum(P)[- (max(s) + 1)])
  robj <- cbind(s = s, "Pr[S>=s]" = pp[s + 1])

  ## give the PMF up to max(s) if verbose = TRUE
  
  if (verbose) {
    attr(robj, "PMF") <- cbind(s = 0L:max(s), "Pr[S=s]" = P)
    if (theta > 0)
      attr(robj, "ELT") <- ELT
  }

  robj
}