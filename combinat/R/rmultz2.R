"rmultz2"<-
function(n, p, draws = length(n))
{
# 19 Feb 1997: From s-news 14 Feb 1997, Alan Zaslavsky
# 11 Mar 1997: Modified by Scott D. Chasalow
#
# Generate random samples from a multinomial(n, p) distn: varying n, 
# fixed p case.
#
	n <- rep(n, length = draws)
	lenp <- length(p)
	tab <- tabulate(sample(lenp, sum(n), TRUE, p) + lenp * rep(1:draws - 1, n),
		nbins = draws * lenp)
	dim(tab) <- c(lenp, draws)
	tab
}

