#' Markov Bound.
#'
#' Function to bound the total losses via the Markov inequality.
#'
#' @param ELT Data frame containing two numeric columns. The column \code{Loss} contains the expected losses from each single occurrence of event. The column \code{Rate} contains the arrival rates of a single occurrence of event.
#' @param s Scalar or numeric vector containing the total losses of interest. 
#' @param t Scalar representing the time period of interest. The default value is \code{t} = 1.
#' @param theta Scalar containing information about the variance of the Gamma distribution: \eqn{sd[X] = x * }\code{theta}. The default value is \code{theta} = 0: the loss associated to an event is considered as a constant.
#' @param cap Scalar representing the financial cap on losses for a single event, i.e. the maximum possible loss caused by a single event. The default value is \code{cap} = Inf.

#' @details Cantelli's inequality states: 
#' \deqn{\Pr( S \geq s) \leq \frac{E[S]}{s}}{Pr( S \ge s) \le E[S]/s,}
#'  @return A numeric matrix, containing the pre-specified losses \code{s} in the first column and the upper bound for the exceedance probabilities in the second column.
#' @keywords markov
#' @export
#' @examples
#' data(UShurricane)
#'
#' # Compress the table to millions of dollars
#'
#' USh.m <- compressELT(ELT(UShurricane), digits = -6)

#' EPC.Markov <- fMarkov(USh.m, s = 1:40)
#' plot(EPC.Markov, type = "l", ylim = c(0, 1))

#' # Assuming the losses follow a Gamma with E[X] = x, and Var[X] = 2 * x

#' EPC.Markov.Gamma <- fMarkov(USh.m, s = 1:40, theta = 2, cap = 5)
#' EPC.Markov.Gamma
#' plot(EPC.Markov.Gamma, type = "l", ylim = c(0, 1))

#' # Compare the two results:

#' plot(EPC.Markov, type = "l", main = "Exceedance Probability Curve", ylim = c(0,1))
#' lines(EPC.Markov.Gamma, col = 2, lty = 2)
#' legend("topright", c("Dirac Delta", expression(paste("Gamma(", 
#' alpha[i] == 1 / theta^2, ", ", beta[i] ==1 / (x[i] * theta^2), ")", " cap =", 5))),  
#' lwd = 2, lty = 1:2, col = 1:2)

fMarkov <- function(ELT, s, t = 1, theta = 0, cap = Inf) {
	
		stopifnot(inherits(ELT, "ELT"),
            s >= 0, 
            is.scalar(t), t > 0,
            is.scalar(theta), theta >= 0,
            is.scalar(cap), cap > 0)

  mom <- getMoments(ELT$Loss, theta = theta, cap = cap, maxmom = 1)
  lambda <- sum(ELT$Rate)
  EY <- (crossprod(ELT$Rate, mom)) / lambda
  
  ES <- lambda * t * EY
  
  pp <- pmin(1, ES / s)
  
  robj <- cbind(s = s, "Upper Pr[S>=s]" = pp)
  
  robj
}