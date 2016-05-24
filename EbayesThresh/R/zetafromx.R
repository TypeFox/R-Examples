"zetafromx" <-
function(xd, cs, pilo=NA, prior = "laplace", a = 0.5)
{
#
#  given a sequence xd, a vector of scale factors cs and
#  a lower limit pilo, find the marginal maximum likelihood
#  estimate of the parameter zeta such that the prior prob
#  is of the form median( pilo, zeta*cs, 1)
#
#  if pilo=NA then it is calculated according to the sample size
#  to corrrespond to the universal threshold
#  
#
#  Find the beta values and the minimum weight if necessary
#
	pr <- substring(prior, 1, 1)
	nx <- length(xd)
	if (is.na(pilo)) pilo <- wfromt(sqrt(2 * log(length(xd))), prior, a)
	if(pr == "l")
		beta <- beta.laplace(xd, a)
	if(pr == "c") beta <- beta.cauchy(xd)
#
#  Find jump points zj in derivative of log likelihood as function
#    of z, and other preliminary calculations
#
	zs1 <- pilo/cs
	zs2 <- 1/cs
	zj <- sort(unique(c(zs1, zs2)))
	cb <- cs * beta
	mz <- length(zj)
	zlmax <- NULL	#
#  Find left and right derivatives at each zj
#   and check which are local minima
#  Check internal zj first
#
	lmin <- rep(FALSE, mz)
	for(j in (2:(mz - 1))) {
		ze <- zj[j]
		cbil <- cb[(ze > zs1) & (ze <= zs2)]
		ld <- sum(cbil/(1 + ze * cbil))
		if(ld <= 0) {
			cbir <- cb[(ze >= zs1) & (ze < zs2)]
			rd <- sum(cbir/(1 + ze * cbir))
			lmin[j] <- (rd >= 0)
		}
	}
#
#  Deal with the two end points in turn, finding right deriv
#   at lower end and left deriv at upper
#
#  In each case the corresponding end point is either a local min
#   or a local max depending on the sign of the relevant deriv
#
	cbir <- cb[zj[1] == zs1]
	rd <- sum(cbir/(1 + zj[1] * cbir))
	if(rd > 0) lmin[1] <- TRUE	else zlmax <- zj[1]
	cbil <- cb[zj[mz] == zs2]
	ld <- sum(cbil/(1 + zj[mz] * cbil))
	if(ld < 0) lmin[mz] <- TRUE else zlmax <- zj[mz]	#
#  Flag all local minima and do a binary search between them to find the local maxima
#
	zlmin <- zj[lmin]
	nlmin <- length(zlmin)
	if(nlmin >= 2) for(j in (2:nlmin)) {
			zlo <- zlmin[j - 1]
			zhi <- zlmin[j]
			ze <- (zlo + zhi)/2
			zstep <- (zhi - zlo)/2
			for(nit in (1:10)) {
				cbi <- cb[(ze >= zs1) & (ze <= zs2)]
				likd <- sum(cbi/(1 + ze * cbi))
				zstep <- zstep/2
				ze <- ze + zstep * sign(likd)
			}
			zlmax <- c(zlmax, ze)
		}
#
#  Evaluate all local maxima and find global max; use smaller value
#   if there is an exact tie for the global maximum.
#
	nlmax <- length(zlmax)
	zm <- rep(NA, nlmax)
	for(j in (1:nlmax)) {
		pz <- pmax(zs1, pmin(zlmax[j], zs2))
		zm[j] <- sum(log(1 + cb * pz))
	}
	zeta <- zlmax[zm == max(zm)]
	zeta <- min(zeta)
	w <- pmin(1, pmax(zeta*cs, pilo) ) 
	return(list(zeta=zeta, w=w, cs=cs, pilo=pilo))
}
