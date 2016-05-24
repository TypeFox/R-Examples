#  SPECIAL FUNCTIONS

logmdigamma <- function(x)
{
#  log(x) - digamma(x)
#  Saves computation of log(x) and avoids subtractive cancellation in digamma(x) when x is large
#  Gordon Smyth, smyth@wehi.edu.au
#  19 Jan 98.  Last revised 9 Dec 2002.
#
	z <- x
	if(any(omit <- is.na(z) | Re(z) <= 0)) {
		ps <- z
		ps[omit] <- NA
		if(any(!omit)) ps[!omit] <- Recall(z[!omit])
		return(ps)
	}
	if(any(small <- Mod(z) < 5)) {
		ps <- z
		x <- z[small]
		ps[small] <- log(x/(x+5)) + Recall(x+5) + 1/x + 1/(x+1) + 1/(x+2) + 1/(x+3) + 1/(x+4)
		if(any(!small)) ps[!small] <- Recall(z[!small])
		return(ps)
	}
	x <- 1/z^2
	tail <- ((x * (-1/12 + ((x * (1/120 + ((x * (-1/252 + ((
		x * (1/240 + ((x * (-1/132 + ((x * (691/32760 + (
		(x * (-1/12 + (3617 * x)/8160)))))))))))))))))))))
	1/(2 * z) - tail
}
