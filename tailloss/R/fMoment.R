#' Moment Bound.
#'
#' Function to bound the total losses via the Moment inequality.
#'
#' @param ELT Data frame containing two numeric columns. The column \code{Loss} contains the expected losses from each single occurrence of event. The column \code{Rate} contains the arrival rates of a single occurrence of event.
#' @param s Scalar or numeric vector containing the total losses of interest. 
#'  @param t Scalar representing the time period of interest. The default value is \code{t} = 1.
#'  @param theta Scalar containing information about the variance of the Gamma distribution: \eqn{sd[X] = x * }\code{theta}. The default value is \code{theta} = 0: the loss associated to an event is considered as a constant.
#' @param cap Scalar representing the financial cap on losses for a single event, i.e. the maximum possible loss caused by a single event. The default value is \code{cap} = Inf.
#'  @param verbose Logical. If \code{TRUE} attaches the minimising index. The default is \code{verbose} = FALSE.

#'  @details Moment inequality states: 
#' \deqn{\Pr(S \geq s) \leq \min_{k = 1, 2 \dots} \frac{E(S^k)}{s^k} }{Pr(S \ge s) \le min_{k = 1, 2 \dots} E(S^k)/s^k }
#' where \eqn{E(S^k)} is the \eqn{k}-th moment of the total loss \eqn{S} distribution.
#'  @return A numeric matrix, containing the pre-specified losses \code{s}  in the first column and the upper bound for the exceedance probabilities in the second column.

#' @keywords moment bound
#' @export
#' @examples
#' data(UShurricane)
#'
#' # Compress the table to millions of dollars
#'
#' USh.m <- compressELT(ELT(UShurricane), digits = -6)

#' EPC.Moment <- fMoment(USh.m, s = 1:40)
#' EPC.Moment
#' plot(EPC.Moment, type = "l", ylim = c(0, 1))

#' # Assuming the losses follow a Gamma with E[X] = x, and Var[X] = 2 * x

#' EPC.Moment.Gamma <- fMoment(USh.m, s = 1:40, theta = 2, cap = 5)
#' EPC.Moment.Gamma
#' plot(EPC.Moment.Gamma, type = "l", ylim = c(0, 1))

#' # Compare the two results:

#' plot(EPC.Moment, type = "l", main = "Exceedance Probability Curve", ylim = c(0, 1))
#' lines(EPC.Moment.Gamma, col = 2, lty = 2)
#' legend("topright", c("Dirac Delta", expression(paste("Gamma(", 
#' alpha[i] == 1 / theta^2, ", ", beta[i] ==1 / (x[i] * theta^2), ")", " cap =", 5))),  
#' lwd = 2, lty = 1:2, col = 1:2)


fMoment <- function(ELT, s, t = 1, theta = 0, cap = Inf, verbose = FALSE) {
	
	stopifnot(inherits(ELT, "ELT"),
            s >= 0, 
            is.scalar(t), t > 0,
            is.scalar(theta), theta >= 0,
            is.scalar(cap), cap > 0,
            is.logical(verbose))
          
    ## Find the highest moment (maxmom) to calculate the moment bound
    ## It is used an algorithm based on the approximation of S with a Gamma distr      
          
    maxmom <- getMaxMom(ELT, theta = theta, cap = cap, s = max(s), t = t)
    
    ## Calculate the first maxmom moments of the true distribution of S 

  	mom <- getMoments(ELT$Loss, theta = theta, cap = cap, maxmom = maxmom)
  	
  	lambda <- sum(ELT$Rate)
  	EY <- drop(crossprod(ELT$Rate, mom)) / lambda
  	
  	ES <- rep(NA, maxmom)
  	lambda <- lambda * t
  	for (k in seq(along = ES)) {
  		ES[k] <- EY[k] # j = 0
  		if (k > 1L) {
  			j <- 1L:(k-1)
  			ES[k] <- ES[k] + sum(choose(k-1, j) * ES[j] * EY[k-j])
  			}
  		ES[k] <- lambda * ES[k]
  	}
  
  ## Calculate the Moment Bound for the true distribution of S 
  ## set k = 1, 2, ... maxmom, and pick the tightest bound
  
  smom <- outer(s, 1:maxmom, "^")
  ppmat <- ES / t(smom)
  pp <- apply(ppmat, 2, min, na.rm = TRUE)
  pp <- pmin(1, pp)
  
  robj <- cbind(s = s, "Upper Pr[S>=s]" = pp)

  ## verbose attaches the minimising index

  if (verbose){
  	khat <- apply(ppmat, 2, which.min)
  	robj <- cbind(robj, khat)
  	}
     
  robj	
  
}