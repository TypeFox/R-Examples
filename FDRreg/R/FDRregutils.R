colnormalize = function(X) {
	# Normalizes each column of x to have mean 0 and sd 1
	Z = X
	d = ncol(X)
	for(j in 1:d) {
		Z[,j] = {X[,j] - mean(X[,j], na.rm=TRUE)}/sd(X[,j], na.rm=TRUE)
	}
	Z;
}

colcenter = function(X) {
	# Normalizes each column of x to have mean 0
	Z = X
	d = ncol(X)
	for(j in 1:d) {
		Z[,j] = {X[,j] - mean(X[,j], na.rm=TRUE)}
	}
	Z;
}

make.bspline.matrix = function(X, order, nknots=5, method="quantile") {
# X = matrix of predictors
# returnval = basis expansion of each column of X in b spline basis
	Z = colnormalize(X)
	d = ncol(X)
	Xs = NULL
	for(j in 1:d) {
		if(method=='quantile') {
			argq = seq(0,1,length=nknots+2)
			argvals = quantile(Z[,j], prob=argq)
		}
		else {
			argvals = seq(min(Z[,j]), max(Z[,j]), length=nknots+2)
		}
		basisj = fda::create.bspline.basis(rangeval=range(argvals), breaks=argvals)
		Xs = cbind(Xs, fda::eval.basis(Z[,j], basisj))
	}
	return(Xs)
}

efron = function(z, nmids=150, pct=-0.01, pct0=0.25, df=10) {
	# estimate f(z) and f_0(z) using Efron (2004)'s method
	N = length(z)
	med = median(z)
	myrange = med + (1 - pct) * (range(z) - med)
	lb = myrange[1]
	ub = myrange[2]
	                
	breaks = seq(lb, ub, length= nmids +1)
	h1 = hist(z, breaks = breaks, plot = FALSE)
	mids = (breaks[-1] + breaks[-(nmids+1)])/2
	zcounts = h1$counts
	glm1 = glm(zcounts ~ splines::ns(mids, df = df), family=poisson)
	zrate = glm1$fit
	D = (zcounts - zrate)/sqrt(zrate+1)
	D = D[-c(1,nmids)]
	if (sum(D^2) > qchisq(0.9, nmids-2-df)) {
		warning(paste0("f(z) misfit = ", round(D, 1), ".  Rerun with increased df."))
	}	
	
	zdens = {zrate/sum(zrate)}/diff(breaks)

	# Now do spline interpolation for the density at the observed points
	ispl2 = splines::interpSpline( zdens ~ mids )
	fz = predict(ispl2, z)$y
	
	# empirical null by central moment matching
	ql = quantile(z, pct0)
	qu = quantile(z, 1-pct0)
	ind0 = intersect(which(z > ql), which(z<qu))
	z0 = z[ind0]
	l0 = log(fz[ind0])
	zmax = z[which.max(l0)]
	lm0 = lm(l0~I(z0-zmax) + I((z0-zmax)^2))
	b0 = coef(lm0)
	sig = as.numeric(sqrt(-1/{2*b0[3]}))
	mu = as.numeric(-b0[2]/(2*b0[3]) + zmax)
	list(mids=mids, breaks=breaks, zcounts=zcounts, zdens=zdens, z=z, fz=fz, mu0=mu, sig0=sig)
}

getFDR = function(postprob) {
	# postprob is a vector of posterior probabilities
	# from which local fdr and (Bayesian) FDR are extracted
	indices = 1:length(postprob)
	iorder = order(postprob, decreasing=TRUE)
	porder = postprob[iorder]
	localfdr.order = 1-porder
	FDR.order = cumsum(localfdr.order)/indices
	localfdr = indices  # placeholder
	localfdr[iorder] = localfdr.order
	FDR = indices  # placeholder
	FDR[iorder] = FDR.order
	
	# Where local fdr is 1, report the most conservative FDR 
	fdrmax = which(localfdr == 1)
	FDR[fdrmax] = max(FDR)
	list(localfdr=localfdr, FDR=FDR)
}

plotFDR = function(fdrr, Q=0.1, showrug=TRUE, showfz=TRUE, showsub=TRUE) {
	N = length(fdrr$z)
	h1 = hist(fdrr$z, fdrr$breaks, plot=FALSE)
	mytitle = paste0('')
	par(mar=c(5,4,1,1))
	mysub = paste0('Grey bars: original z scores\nRed bars: fraction signals in each bin')
	plot(h1, col='lightgrey', border='grey', main=mytitle, xlab='', ylab='', ylim=range(h1$counts), axes=FALSE)
	axis(1, pos=0, tick=FALSE)
	axis(2, tick=FALSE, las=1)
	zcut = data.frame(prob=fdrr$postprob, bucket=cut(fdrr$z, fdrr$breaks))
	pmean = mosaic::maggregate(prob~bucket, data=zcut, FUN=mean)
	pmean[is.na(pmean)] = 0
	par(new=TRUE)
	h2 = h1
	h2$counts = h2$counts * pmean
	plot(h2, col='red', border='grey', ylim=range(h1$counts), axes=FALSE, xlab='', ylab='', main='')
	lines(fdrr$grid, (fdrr$grid.f0z/sum(fdrr$grid.f0z))*N, col='blue', lty='solid', lwd=1)
	if(showfz) {
		lines(fdrr$grid, (fdrr$grid.fz/sum(fdrr$grid.fz))*N, col='black')
		legend('topright', c(expression(f(z)), expression(f[0](z))), lty=c('solid', 'solid'), col=c('black', 'blue'), bty='n')
	}
	else {
		legend('topright', c(expression(f[0](z))), lty=c('solid'), col=c('blue'), bty='n')	
	}
	if(showrug) {
		rug( fdrr$z[fdrr$FDR < Q], ticksize=0.03, col='black')
		mysub = paste0(mysub, '\nBlack rug: discoveries at FDR = ', Q)
	}
	if(showsub) title(sub=mysub)
	par(new=FALSE)
}
