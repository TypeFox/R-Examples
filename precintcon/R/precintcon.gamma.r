#' @noRd
#' @name precintcon.gamma
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.gamma 
#' @title Gamma distribution fit
#' @description Fits a precipitation serie to a gamma distribution. 
#' @usage precintcon.gamma(p, period) 
#' @param p is a cumulative precipition serie aggregated by a \code{period}.
#' @param period is the aggregation period.
#' @return A fitted gamma distribution. 
#' @seealso
#' \code{\link{spi}}
#' \code{\link{read.data}}
#' \code{\link{as.daily}}
#' \code{\link{as.monthly}} 
#' @keywords precipitation stardardized precipitation index 
precintcon.gamma <- function(p, period) {
	
	if (length(p) <= 0)
		stop("Invalid input data. The size of p should be greater than 0.");
	
	gamma_ <- data.frame()

	for (j in 1:12) {

		im = (period + j - 1) %% 12 + 1
		
		q <- p[seq(j, length(p) + period, by=12)]
		q <- q[!is.na(q)]

		pzero <- sum(q==0) / length(q)
		
		avg <- mean(q[q > 0.0])

		alpha <- 0.0
		beta  <- avg
		gamm  <- 1.0
		
		pgz <- length(q[q > 0.0])
		
		if ( pgz >= 1) {
			alpha <- log(avg) - sum(log(q[q > 0.0])) / pgz 
			gamm <- (1.0 + sqrt(1.0 + 4.0 * alpha / 3.0)) / (4.0 * alpha)
			beta  <- avg / gamm
		} 
				
		gamma_ <- rbind(gamma_, data.frame(month=im, alpha=alpha, beta=beta, gamm=gamm, pzero=pzero))
	}
	
	o <- order(gamma_$month)
	
	gamma_ <- gamma_[o, ]
	
	return(gamma_)
}