"isotone" <-
function(x, wt = rep(1, length(x)), increasing = FALSE)
{
#
#   find the weighted least squares isotone fit to the 
#   sequence x, the weights given by the sequence wt
#
#   if increasing == TRUE the curve is set to be increasing, 
#   otherwise to be decreasing
#
#   the vector ip contains the indices on the original scale of the
#   breaks in the regression at each stage
#
	nn <- length(x)
	if(nn == 1)
		return(x)
	if(!increasing)
		x <-  - x
	ip <- (1:nn)
	dx <- diff(x)
	nx <- length(x)
	while((nx > 1) && (min(dx) < 0)) {
#
#  do single pool-adjacent-violators step
#
#  find all local minima and maxima
#
		jmax <- (1:nx)[c(dx <= 0, FALSE) & c(TRUE, dx > 0)]
		jmin <- (1:nx)[c(dx > 0, TRUE) & c(FALSE, dx <= 0)]	#	
#  do pav step for each pair of maxima and minima
#
#  add up weights within subsequence that is pooled
#  set first element of subsequence to the weighted average
#  the first weight to the sum of the weights within the subsequence
#    and remainder of the subsequence to NA
#
		for(jb in (1:length(jmax))) {
			ind <- (jmax[jb]:jmin[jb])
			wtn <- sum(wt[ind])
			x[jmax[jb]] <- sum(wt[ind] * x[ind])/wtn
			wt[jmax[jb]] <- wtn
			x[(jmax[jb] + 1):jmin[jb]] <- NA
		}
#
#  clean up within iteration, eliminating the parts of sequences that were set
#   to NA
#
		ind <- !is.na(x)
		x <- x[ind]
		wt <- wt[ind]
		ip <- ip[ind]
		dx <- diff(x)
		nx <- length(x)
	}
# 
#  final cleanup:  reconstruct z at all points by repeating the pooled values
#    the appropriate number of times
#
	jj <- rep(0, nn)
	jj[ip] <- 1
	z <- x[cumsum(jj)]
	if(!increasing)
		z <-  - z
	return(z)
}
