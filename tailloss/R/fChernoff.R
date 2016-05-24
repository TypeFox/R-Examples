#' Chernoff Bound.
#'
#' Function to bound the total losses via the Chernoff inequality.
#'
#' @param ELT Data frame containing two numeric columns. The column \code{Loss} contains the expected losses from each single occurrence of event. The column \code{Rate} contains the arrival rates of a single occurrence of event.
#' @param s Scalar or numeric vector containing the total losses of interest. 
#' @param t Scalar representing the time period of interest. The default value is \code{t} = 1.
#' @param theta Scalar containing information about the variance of the Gamma distribution: \eqn{sd[X] = x * }\code{theta}. The default value is \code{theta} = 0: the loss associated to an event is considered as a constant.
#' @param cap Scalar representing the financial cap on losses for a single event, i.e. the maximum possible loss caused by a single event. The default value is \code{cap} = Inf.
#' @param nk Number of optimisation points.
#' @param verbose Logical. If \code{TRUE} attaches the minimising index. The default is \code{verbose} = FALSE.

#' @details Chernoff's inequality states: 
#' \deqn{\Pr(S \geq s) \leq \inf_{k > 0} e^{-k s} M_S(k) }{Pr(S \ge s) \le inf_{k > 0} e^{-k s} M_S(k) }
#' where \eqn{M_S(k)} is the Moment Generating Function (MGF) of the total loss S.
#' The \code{fChernoff} function optimises the bound over a fixed set of \code{nk} discrete values.
#'  @return A numeric matrix, containing the pre-specified losses \code{s} in the first column and the upper bound for the exceedance probabilities in the second column.

#' @keywords chernoff
#' @export
#' @examples
#' data(UShurricane)
#'
#' # Compress the table to millions of dollars
#'
#' USh.m <- compressELT(ELT(UShurricane), digits = -6)

#' EPC.Chernoff <- fChernoff(USh.m, s = 1:40)
#' EPC.Chernoff
#' plot(EPC.Chernoff, type = "l", ylim = c(0, 1))

#' # Assuming the losses follow a Gamma with E[X] = x, and Var[X] = 2 * x

#' EPC.Chernoff.Gamma <- fChernoff(USh.m, s = 1:40, theta = 2, cap = 5)
#' EPC.Chernoff.Gamma
#' plot(EPC.Chernoff.Gamma, type = "l", ylim = c(0, 1))

#' # Compare the two results:

#' plot(EPC.Chernoff, type = "l", main = "Exceedance Probability Curve", ylim = c(0, 1))
#' lines(EPC.Chernoff.Gamma, col = 2, lty = 2)
#' legend("topright", c("Dirac Delta", expression(paste("Gamma(", 
#' alpha[i] == 1 / theta^2, ", ", beta[i] ==1 / (x[i] * theta^2), ")", " cap =", 5))),  
#' lwd = 2, lty = 1:2, col = 1:2)

fChernoff <- function(ELT, s, t = 1, theta = 0, cap = Inf, nk = 1001, verbose = FALSE){
	
		stopifnot(inherits(ELT, "ELT"),
            s >= 0, 
            is.scalar(t), t > 0,
            is.scalar(theta), theta >= 0, 
            is.scalar(cap), cap > 0,
            is.scalar(nk), nk > 0, nk == floor(nk),
            is.logical(verbose))
  
    lambda <- sum(ELT$Rate)
    EY <- sum(ELT$Rate * ELT$Loss) / lambda
    
    ES <- lambda * t * EY

    pp <- rep(1, length(s))
    stillin <- s > ES
  
    if (theta == 0) {
    	
    	if (!is.infinite(cap))
    			ELT$Loss <- pmin(cap, ELT$Loss)
    			
    	k <- seq(from = 0, to = min(100 / ELT$Loss), length = nk)
      
      ## going to minimise the log of the optimand

      ppmat <- lambda * t * (colSums(ELT$Rate * exp(ELT$Loss %o% k)) / lambda - 1)
      ppmat <- (-k) %o% s[stillin] + drop(ppmat)
      ppmat <- exp(ppmat)
      pp[stillin] <- apply(ppmat, 2, min)
      
      } else { # theta > 0
      	
      alpha <- 1 / theta^2
      beta <- alpha / ELT$Loss
      	  
      k <- seq(from = 0, to = min(100 / ELT$Loss, beta), length = nk + 1L)[- (nk + 1L)]
      	 
      ## going to minimise the log of the optimand

	if (is.infinite(cap)){ # 
		
		ppmat <- lambda * t * (
		colSums(outer(ELT$Rate, k, function(h1,h2){
    	h1 * (beta / (beta - h2))^alpha
    	})) / lambda - 1)
    	
    	ppmat <- (-k) %o% s[stillin] + drop(ppmat)
    		ppmat <- exp(ppmat)
    		pp[stillin] <- apply(ppmat, 2, min)
		
	} else { # cap < Inf

		u <- cap
		p <- pgamma(u, shape = alpha, rate = beta)
	
		ppmat <- lambda * t * ((
			colSums(outer(ELT$Rate, k, function(h1, h2){
    		h1 * (beta / (beta - h2))^alpha * igamma(alpha, (beta - h2) * u) / gamma(alpha)
    		})) + 
    		sum((1 - p) * ELT$Rate) * exp(u * k)
    		) / lambda - 1 )
    
    	ppmat <- (-k) %o% s[stillin] + drop(ppmat)
    	ppmat <- exp(ppmat)
    	pp[stillin] <- apply(ppmat, 2, min)
      }
    
    }
    
   robj <- cbind(s = s, "Upper Pr[S>=s]" = pp)
  
  ## verbose attaches the minimising index

  if (verbose) {
      khat <- rep(0, length(s))
      khat[stillin] <- k[apply(ppmat, 2, which.min)]
      robj <- cbind(robj, khat = khat)
  }

  robj
}