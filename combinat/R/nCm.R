"nCm"<-
function(n, m, tol = 9.9999999999999984e-009)
{
#  DATE WRITTEN:  7 June 1995               LAST REVISED:  10 July 1995
#  AUTHOR:  Scott Chasalow
#
#  DESCRIPTION: 
#        Compute the binomial coefficient ("n choose m"),  where n is any 
#        real number and m is any integer.  Arguments n and m may be vectors;
#        they will be replicated as necessary to have the same length.
#
#        Argument tol controls rounding of results to integers.  If the
#        difference between a value and its nearest integer is less than tol,  
#        the value returned will be rounded to its nearest integer.  To turn
#        off rounding, use tol = 0.  Values of tol greater than the default
#        should be used only with great caution, unless you are certain only
#        integer values should be returned.
#
#  REFERENCE: 
#        Feller (1968) An Introduction to Probability Theory and Its 
#        Applications, Volume I, 3rd Edition, pp 50, 63.
#
	len <- max(length(n), length(m))
	out <- numeric(len)
	n <- rep(n, length = len)
	m <- rep(m, length = len)
	mint <- (trunc(m) == m)
	out[!mint] <- NA
	out[m == 0] <- 1	# out[mint & (m < 0 | (m > 0 & n == 0))] <-  0
	whichm <- (mint & m > 0)
	whichn <- (n < 0)
	which <- (whichm & whichn)
	if(any(which)) {
		nnow <- n[which]
		mnow <- m[which]
		out[which] <- ((-1)^mnow) * Recall(mnow - nnow - 1, mnow)
	}
	whichn <- (n > 0)
	nint <- (trunc(n) == n)
	which <- (whichm & whichn & !nint & n < m)
	if(any(which)) {
		nnow <- n[which]
		mnow <- m[which]
		foo <- function(j, nn, mm)
		{
			n <- nn[j]
			m <- mm[j]
			iseq <- seq(n - m + 1, n)
			negs <- sum(iseq < 0)
			((-1)^negs) * exp(sum(log(abs(iseq))) - lgamma(m + 1))
		}
		out[which] <- unlist(lapply(seq(along = nnow), foo, nn = nnow, 
			mm = mnow))
	}
	which <- (whichm & whichn & n >= m)
	nnow <- n[which]
	mnow <- m[which]
	out[which] <- exp(lgamma(nnow + 1) - lgamma(mnow + 1) - lgamma(nnow - 
		mnow + 1))
	nna <- !is.na(out)
	outnow <- out[nna]
	rout <- round(outnow)
	smalldif <- abs(rout - outnow) < tol
	outnow[smalldif] <- rout[smalldif]
	out[nna] <- outnow
	out
}

