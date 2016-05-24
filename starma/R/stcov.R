# The 'stcov' function is defined as, per (Pfeifer & Stuart, 1980), the 
# following expression:
# 	cov(l,k,s) = E[ (W(l)*z(t))' * (W(k)*z(t+s)) ] / N
#		     = Tr( W(k)'*W(l) * E[z(t)*z(t+s)'] ) / N

stcov <- function(data, wlist, slag1, slag2, tlag) {

	if (is.matrix(wlist))
		wlist <- list(diag(dim(wlist)[1]), wlist)
	
	if (is.data.frame(data))
		return(stcovCPP(as.matrix(data), wlist, slag1, slag2, tlag))
	else
		return(stcovCPP(data, wlist, slag1, slag2, tlag))

}

# Ideas to optimize code:
# - Try to export data.frame into C++, since converting it into a matrix
#   beforehand takes a lot of time.
# - Try to express wlist in a RcppArmadillo native code. As it is, it must
#   create two matrices w1 and w2.
# - Using sparse matrices for weight matrices.