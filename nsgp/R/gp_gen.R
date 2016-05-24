



simgp_clust = function(n.clusters=30, n.outliers=50, t=c(0,10) ) {
	ts = seq(t[1],t[2], (t[2]-t[1])/100)
	n = length(ts)
	n.gp = n.clusters + n.outliers
	
	K.noise = .RBF(ts,ts,0.3,5)
	
	# draw lengthscales and sigma_f's
	ells = rgamma(n.gp,4,scale=1)
	sfs = rgamma(n.gp,4,scale=0.25)
	reps = round(runif(30,5,10))
	
	obspoints = seq(1,n,8)
	truegps = list()
	data = mat.or.vec(sum(reps),length(obspoints))
	ind = 1
	
	# draw 30 clusters with 5-10 reps
	for (i in 1:n.clusters) {
		f.mean = mvrnorm(1, rep(0,n), .RBF(ts,ts,sfs[i],ells[i])) # draw sample
		gp = gpsimple(ts, f.mean, K.noise, rep(0,n)) # add noise
		
		ys = mvrnorm(reps[i], gp$mean, gp$cov)
		
		truegps[[i]] = gp
		
		data[ind:(ind+reps[i]-1),] = ys[,obspoints]
		ind = ind + reps[i]
#		plot(gp); for (j in 1:rcount) { lines(ts, reps[j,]) }
	}
	
	plot(NA, xlim=c(0,10), ylim=c(-5,5))
	for (i in 1:sum(reps)) {
		lines(ts[obspoints], data[i,])
	}

	# draw 50 outlier curves
	for (i in 1:n.outliers) {
		f.mean = mvrnorm(1, rep(0,n), .RBF(ts,ts,sfs[i],ells[i])) # draw sample
		gp = gpsimple(ts, f.mean, K.noise, rep(0,n)) # add noise
		rep = mvrnorm(1, gp$mean, gp$cov)
		
		truegps[[i+n.clusters]] = gp
	}
}

#' @export
#' @title Generate simulated GP models
#' @param N number of gp's 
#' @param reps replicate observations 
#' @param obs observation timepoints 
#' @param xs target timepoints 
#' @param filename file to save the results 
#' @param l.noise Noise model lengthscale 
#' @param l.shape shape of lengthscale Gamma distribution 
#' @param l.scale scale of lengthscale Gamma distribution 
#' @param sigmaf.noise Noise model sigma.f 
#' @param sigmaf.shape shape of sigma.f Gamma distribution 
#' @param sigmaf.scale scale of sigma.f Gamma distribution 
#' @return List with 
#'  \item{simdata}{Simulated datamatrix}
#'  \item{gps}{Simulated GPs}
#' @seealso \code{\link{simulategp.perturbed}}
simulategp = function(N=100, reps=3, obs=c(0,5,10,15,20), xs=seq(0,20,0.2), filename=NULL,
								l.noise=12, l.shape=4, l.scale=4,
								sigmaf.noise=0.25, sigmaf.shape=1.5, sigmaf.scale=0.3) {
	n = length(xs)
	K.noise = .RBF(xs,xs,sigmaf.noise,l.noise)	
	ells = rgamma(N,l.shape,scale=l.scale)
	sigmafs = rgamma(N,sigmaf.scale,scale=sigmaf.scale)
	
	gps = list()
	for (i in 1:N) {
		f.mean = mvrnorm(1, rep(0, n), .RBF(xs,xs,sigmafs[i],ells[i])) # draw sample
		gp = gpsimple(xs, f.mean, K.noise, rep(0, n))
		gps[[i]] = gp
#		plot(gp, samples=3)
	}
	
	simdata = data.frame(mat.or.vec(N, length(obs)*reps))
	names(simdata) = rep(obs,reps)
	for (i in 1:N) {
		res = mvrnorm(reps, gps[[i]]$mean, gps[[i]]$cov)
		simdata[i,] = as.vector(t(res[,match(obs,xs)]))
	}

	if (!is.null(filename)) {
		save(simdata,gps, file=filename)
	}
	
	return(list("simdata"=simdata,"gps"=gps))			
}

plotdistr = function(l.shape=4, l.scale=4, sigmaf.shape=1.5, sigmaf.scale=0.3) {
	maxs = sum(dgamma(0:1000,sigmaf.shape,scale=sigmaf.scale) > 0.001)
	maxl = sum(dgamma(0:1000,l.shape,scale=l.scale) > 0.001)
	
	sfs = seq(0,maxs,0.01)
	ls = seq(0,maxl,0.01)

	par(mfrow=c(1,2))
	plot(sfs, dgamma(sfs,sigmaf.shape,scale=sigmaf.scale), 
		  t='l', ylab='Density', xlab='sigma_f', main='sigma_f distribution')
	plot(ls, dgamma(ls,l.shape,scale=l.scale), 
		  t='l', ylab='Density', xlab='ell', main='ell distribution')
}


#' @export
#' @title Generate simulated perturbed GP models
#' @description Takes pregenerated simulated GP models as input, and
#'  models both the perturbation and additional perturbation variance as
#'   GP models, that are all combined into a perturbed GP model. We sample 
#'   \code{reps} data points from this at timepoints \code{obs}
#' @param N number of gp's
#' @param reps replicate observations 
#' @param obs observation timepoints 
#' @param xs target timepoints 
#' @param filename file to save the results 
#' @param l.noise Noise model lengthscale 
#' @param l.shape shape of lengthscale Gamma distribution 
#' @param l.scale scale of lengthscale Gamma distribution 
#' @param sigmaf.noise Noise model sigma.f 
#' @param sigmaf.shape shape of sigma.f Gamma distribution 
#' @param sigmaf.scale scale of sigma.f Gamma distribution 
#' @param gps.ctrl the GPS models to be perturbed
#' @return List with
#'  \item{simdata}{Simulated datamatrix}
#'  \item{gps}{Simulated GPs}
#' @seealso \code{\link{simulategp}}
simulategp.perturbed = function(N=100, reps=3, obs=c(0,5,10,15,20), 
										  gps.ctrl=NULL, filename=NULL, xs=seq(0,20,0.2),
										  l.noise=12, l.shape=2, l.scale=2.5,
										  sigmaf.noise=0.25, sigmaf.shape=4, sigmaf.scale=0.5) {
	n = length(xs)
	zeros = rep(0, n)
	ones = rep(1,n)
	K.noise = .RBF(xs,xs,sigmaf.noise,l.noise)	
	
	# case distributions
	draw.lmax = function(n=1) { rgamma(n,l.shape,scale=l.scale) }
	draw.lmin = function(n=1,lmax) { runif(n, min(0.5,lmax), lmax) }
	draw.sigmaf = function(n=1) { rgamma(n,sigmaf.shape,scale=sigmaf.scale) }
	draw.a = function(n=1) { runif(n, min(xs),max(xs)/2) }
	draw.b = function(n=1,a,ell) { runif(n, a+(ell*2), max(max(xs),a+(ell*2))) }
	draw.c = function(n=1) { exp(runif(n, -2, 3)) }
	
	ells.case = draw.lmax(N)
	ells.case[ells.case<1] = 1
	ells.min.case = draw.lmin(N, ells.case)
	ells.min.case[ells.min.case<1] = 1
	sigmafs.case = draw.sigmaf(N)
	sigmafs.case[sigmafs.case<0.5] = 0.5
	cs.case = draw.c(N)
	
	gps.case = list()
	intervals = list()
		
	for (i in 1:N) {
		# generate differential time periods
		inds = rep(FALSE,n)
		points = c()
		a = draw.a()
		b = draw.b(a=a, ell= .bexp(a, ells.case[i], ells.min.case[i], cs.case[i]) )
		inds = inds | (xs>=a & xs<=b)
		intervals[[i]] = c(a,b)
		
		# control GP
		gpctrl = gps.ctrl[[i]]
		
		# control curve from the previously generated "steady-state" GP
		f.mean.ctrl = gps.ctrl[[i]]$mean
		
		# a perturbation distribution
		gpcase.prior = gpr_posterior(xs[!inds], zeros[!inds],xs,rep(0.01,n), 
											  function(x1,x2) .RBF(x1,x2,sigmafs.case[i],ells.case[i]))
		# draw a curve from it
		f.mean.case = mvrnorm(1, gpcase.prior$mean, gpcase.prior$cov)
		
		# additional noise GP
		gpnoise = gpr_posterior(xs[!inds], zeros[!inds],xs,rep(0.01,n), function(x1,x2) .RBF(x1,x2,0.1,2))
		
		# final perturbed GP
		gpcase = gpsimple(xs, f.mean.ctrl + f.mean.case, gps.ctrl[[i]]$cov + gpnoise$cov, zeros)
		gps.case[[i]] = gpcase

		# nice example plot
#		pdf('~/Work/Rosiris/figures/simdata-example.pdf', 11,3.5)
#		par(mfrow=c(1,4))
#		plot(gpctrl, title='Control GP')
#		plot(gpcase.prior, title='Differential interval posterior + sample', plotdata=F)
#		rect(a, -10, b,10, col=adjustcolor('grey',0.10))
#		lines(xs, f.mean.case)
#		plot(gpnoise, title='Additional case variance', plotdata=F)
#		plot(gpcase, title='Perturbed GP')
	}
	
	# differential case data
	casedata = mat.or.vec(N, length(obs)*reps)
	for (i in 1:N) {
		res = mvrnorm(reps, gps.case[[i]]$mean, gps.case[[i]]$cov)
		casedata[i,] = as.vector(t(res[,match(obs,xs)]))
	}

	if (!is.null(filename)) {
		save(casedata,gps.case,intervals, file=filename)
	}	
	
	return(list("casedata"=casedata,"casegps"=gps.case,"caseintervals"=intervals))
}
