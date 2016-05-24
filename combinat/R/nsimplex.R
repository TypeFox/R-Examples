"nsimplex"<-
function(p, n)
{
# DATE WRITTEN:  24 Dec 1997 		 LAST REVISED:  24 Dec 1997
# AUTHOR:  Scott D. Chasalow  (Scott.Chasalow@users.pv.wau.nl)
#
# DESCRIPTION:
#       Computes the number of points on a {p, n}-simplex lattice; that is, the
#	number of p-part compositions of n. This gives the number of points in
#	the support space of a Multinomial(n, q) distribution, where
#	p == length(q).
#
#	Arguments p and n are replicated as necessary to have the length of the
#	longer of them.
#
# REQUIRED ARGUMENTS:
#	p	vector of (usually non-negative) integers
#	n	vector of (usually non-negative) integers
# 
	mlen <- max(length(p), length(n))
	p <- rep(p, length = mlen)
	n <- rep(n, length = mlen)
	out <- nCm(n + p - 1, n)
	out[p < 0] <- 0
	out
}

