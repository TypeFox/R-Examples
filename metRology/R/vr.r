#===============================================
# ONE-WAY RANDOM EFFECTS MODEL MLE
# 1. USING SCORE ESTIMATING EQUATIONS
# 2. USING THE VANGEL-RUKHIN ITERATIVE ALGORITHM
#
# NIST
# Hugo Gasca-Aragon
# created: May 2009
# last update:
#	Apr 2010
#	The VR algorithm was updated with
#     a Gauss-Newton search for the roots
#     instead of the linear seach
#
#	Jun 2010
#	The VR result was extended to include
#	if the convergence criteria was met,
#	if the reduced model is suggested
#	(no evidence of random effect)
#===============================================

#Changes:
#
# 2014-03-04  Removed "require(MASS)" (deprecated) - SLRE
#

#
# Suggestion [SLRE]
# i) Add construct.loc.estimate and return location estimate
#

mle.1wre<- function(x, s2, n, init.mu=mean(x), init.sigma2=var(x), labels=c(1:length(x)), max.iter=200, tol=.Machine$double.eps^0.5, trace=FALSE)
{
#
# parameters
#
# x=reported sample means, s2=reported sample variances, n=reported sample sizes, labels=participant labels
#

	# sort the datapoints by labels
	labels <- as.character(labels)	#Added SLRE 2010-08-01
	x<- x[order(labels)]
	s2<- s2[order(labels)]
	n<- n[order(labels)]
	labels<- labels[order(labels)]

	# remove all datapoints with undefined s2
	x<- x[!is.na(s2)]
	n<- n[!is.na(s2)]
	labels<- labels[!is.na(s2)]
	s2<- s2[!is.na(s2)]

	p<- length(x)
	N<- sum(n)
	Theta<- matrix(0, max.iter+1, p+2)

	# set the initial values for sigma2 and sigmai2
	mu<- init.mu
	sigma2<- init.sigma2
	sigmai2<- s2
	t<- 1

	Theta[t,]<- c(mu, sigmai2, sigma2)
	cur.rel.abs.error<- Inf
	llh<- 0

	while ((cur.rel.abs.error> tol) && (t<max.iter) && all(sigmai2>0) && (sigma2>0)) {

		mu<- Theta[t,1]
		sigmai2<- Theta[t,c(2:(p+1))]
		sigma2<- Theta[t,p+2]
		wi<- n/(sigmai2+n*sigma2)
		var.mu<- 1/sum(wi)

		llh<- -N/2*log(2*pi)-sum((n-1)*log(sigmai2))/2-sum(log(sigmai2+n*sigma2))/2-sum((n-1)*s2/sigmai2)/2-sum(n*(x-mu)^2/(sigmai2+n*sigma2))/2
		if (trace) {
			print( c(t=t, mu=mu, var.mu=var.mu, sigma2=sigma2, sigmai2=sigmai2, llh=llh) )
		}

		# get the score vector
		S<- c(sum((x-mu)/(sigma2+sigmai2/n)), -(n-1)/sigmai2/2-1/(n*sigma2+sigmai2)/2+(n-1)*s2/sigmai2^2/2+n*(x-mu)^2/(n*sigma2+sigmai2)^2/2, -sum(1/(sigma2+sigmai2/n))/2+sum((x-mu)^2/(sigma2+sigmai2/n)^2)/2)

		if (trace) {
			print( "Score vector:" )
			print( S )
		}

		# get the information matrix
		I<- matrix(0, p+2, p+2)

		I[1,1]<- sum(1/(sigma2+sigmai2/n))
		for (j in 1:p) {
			I[j+1,j+1]<- (1/2)*(n[j]-1)/sigmai2[j]^2+1/(sigmai2[j]+n[j]*sigma2)^2/2
			I[j+1,p+2]<- n[j]/2/(sigmai2[j]+n[j]*sigma2)^2
			I[p+2,j+1]<- I[j+1,p+2]
		}
		I[p+2,p+2]<- sum(1/(sigma2+sigmai2/n)^2)/2


		# get the inverse of the information matrix
		Iinv<- ginv(I)

		# update the estimates
		Theta[t+1,]<- Theta[t,] + Iinv %*% S

		# if the new estimate of sigma is negative set to zero
		#Theta[t+1,c(FALSE,Theta[t+1,2:(p+2)]<0)]<- 0

		# compute the maximum relative absolute error
		cur.rel.abs.error<- max(abs((Theta[t+1,]-Theta[t,])/Theta[t,]))
		if (trace) print( c("current relative absolute error <=", cur.rel.abs.error) )

		t<- t+1

		mu<- Theta[t,1]
		sigmai2<- Theta[t,c(2:(p+1))]
		sigma2<- Theta[t,p+2]

		if (is.na(sigma2)) stop("sigma2 became undefined.")
		if (any(is.na(sigmai2))) stop("some sigmai2 became undefined")
		if (is.na(mu)) stop("mu became undefined")
		if (is.na(cur.rel.abs.error)) stop("current relative absolute error became undefined")
	} # while

	if ((t==max.iter)&&(cur.rel.abs.error>tol)) {
		warning("Non convergence or slow convergence condition was found.")
	} 

	mu<- Theta[t,1]
	sigmai2<- Theta[t,c(2:(p+1))]
	sigma2<- Theta[t,p+2]
	wi<- n/(sigmai2+n*sigma2)
	var.mu<- 1/sum(wi)
	llh<- -N/2*log(2*pi)-sum((n-1)*log(sigmai2))/2-sum(log(sigmai2+n*sigma2))/2-sum((n-1)*s2/sigmai2)/2-sum(n*(x-mu)^2/(sigmai2+n*sigma2))/2
	if (trace) print( c(t=t, mu=mu, var.mu=var.mu, sigma2=sigma2, sigmai2=sigmai2, llh=llh) )

	result<- list( mu=as.vector(mu), var.mu=as.vector(var.mu), sigma2=as.vector(sigma2), llh=as.vector(llh), tot.iter=as.vector(t), max.rel.abs.error=as.vector(cur.rel.abs.error), sigmai2=as.vector(sigmai2) )
	class(result)<- "summary.mle.1wre"

	return(result)

} # function block



.newton.raphson<- function(f, fp, init.value, max.tol=.Machine$double.eps^0.5, trace=FALSE, range.tol=0) {
	max.iter<- 40
	tol<- 1
	iter<- 0
	x<- init.value
	while ((iter < max.iter)&(tol>max.tol)) {
		xn<- x-f(x)/fp(x)
		if (trace) print( c(iter, xn, x, f(x), fp(x), tol) )
		iter<- iter+1

		if (range.tol==0) {
			tol<- abs(xn-x)
		} else 
		if (range.tol==1) {
			tol<- abs(xn-x)/abs(max(x,xn))
		} else {
			tol<- abs(f(xn))
		}

		x<- xn
		if (trace) print( c(iter, xn, x, f(x), fp(x), tol) )
	}
	if (x>1) x<- 1
	if (x<0) x<- 0
	if (trace) print( c(iter, xn, x, f(x), fp(x), tol) )
	return(x)
}




vr.mle<- function(x, s2, n, init.mu=mean(x), init.sigma2=var(x), labels=c(1:length(x)), max.iter=1000, tol=.Machine$double.eps^0.5, trace=FALSE)
{
#
# parameters
#
# x=reported sample means, s2=reported sample variances, n=reported sample sizes, labels=participant labels
#

	# sort the datapoints by labels
	labels <- as.character(labels)	#Added SLRE 2010-08-01
	x<- x[order(labels)]
	s2<- s2[order(labels)]
	n<- n[order(labels)]
	labels<- labels[order(labels)]

	# remove all datapoints with undefined s2
	x<- x[!is.na(s2)]
	n<- n[!is.na(s2)]
	labels<- labels[!is.na(s2)]
	s2<- s2[!is.na(s2)]

	p<- length(x)
	N<- sum(n)
	Theta<- matrix(0, max.iter+1, p+2)

	# set the initial values for sigma2, sigmai2 and gammai
	mu<- init.mu
	sigma2<- init.sigma2
	sigmai2<- s2

	gammai<- sigma2/(sigma2+sigmai2/n)

	t<- 1

	Theta[t,]<- c(mu, gammai, sigma2)
	cur.rel.abs.error<- Inf
	llh<- 0

	while ((cur.rel.abs.error> tol) && (t<max.iter) && all(gammai>0) && (sigma2>0)) {

		mu<- Theta[t,1]
		gammai<- Theta[t,c(2:(p+1))]
		sigma2<- Theta[t,p+2]
		var.mu<- sigma2/sum(gammai)

		llh<- -N/2*log(2*pi)-sum(n*log(n))/2+sum(n*log(gammai/sigma2))/2-sum((n-1)*log(1-gammai))/2-sum(gammai*((x-mu)^2+(n-1)*s2/n/(1-gammai)))/sigma2/2

		# Iterative update
		ai<- sigma2/(x-mu)^2
		bi<- s2/n/(x-mu)^2

		# solve for weights gammai
		for (i in 1:p) {
			#print( c(i=i, p=Theta[t, 1+i], steptol=tol, gammai[i], ai[i], bi[i], n[i]) )

			nlm.res<- .newton.raphson(function(x) {x^3-(ai[i]+2)*x^2+((n[i]+1)*ai[i]+(n[i]-1)*bi[i]+1)*x-n[i]*ai[i]}, function(x) {3*x^2-2*(ai[i]+2)*x+(n[i]+1)*ai[i]+(n[i]-1)*bi[i]+1}, init.value=Theta[t,1+i], max.tol=tol, trace=trace)
			# nlm(f=function(x) {x^4/4-(ai[i]+2)*x^3/3+((n[i]+1)*ai[i]+(n[i]-1)*bi[i]+1)*x^2/2-n[i]*ai[i]*x}, p=c(Theta[t,1+i]), steptol=tol)
			# print( c(code=nlm.res$code) )
			gammai[i]<- nlm.res #$estimate
			if (gammai[i]>1) stop("gamma[",i,"]>1")
			if (gammai[i]<0) stop("gamma[",i,"]<0")
		}
		
		mu<- sum(gammai*x)/sum(gammai)
		sigma2<- sum(gammai^2*(x-mu)^2)/sum(gammai)

		Theta[t+1,]<- c(mu, gammai, sigma2)

		# compute the maximum relative absolute error
		cur.rel.abs.error<- max(abs((Theta[t+1,]-Theta[t,])/Theta[t,]))

		t<- t+1

		mu<- Theta[t,1]
		gammai<- Theta[t,c(2:(p+1))]
		sigma2<- Theta[t,p+2]

		if (is.na(sigma2)) stop("sigma2 became undefined.")
		if (any(is.na(gammai))) stop("some gammai became undefined")
		if (is.na(mu)) stop("mu became undefined")
		if (is.na(cur.rel.abs.error)) stop("current relative absolute error became undefined")
		if (trace) print( c(t, Theta[t,], llh) )
	} # while

	converged<- TRUE
	if ((t==max.iter)&&(cur.rel.abs.error>tol)) {
		# convergence condition is met
		converged<- FALSE
	} 

	mu<- Theta[t,1]
	gammai<- Theta[t,c(2:(p+1))]
	sigma2<- Theta[t,p+2]
	var.mu<- sigma2/sum(gammai)

	# iterative reweighted mean
	if (sigma2==0) {
		si<- sqrt(s2)
		sigmai<- si
		iter<- 0
		tol<- 1
		mu<- mean(x)
		max.tol<- .Machine$double.eps^0.5
		while ((iter<100) & (tol>max.tol)) {
			mu.n<- sum(x/sigmai^2)/sum(1/sigmai^2)
			sigmai.n<- sqrt(n*(x-mu)^2+(n-1)*si^2)/n
			tol<- abs(mu.n-mu)
			mu<- mu.n
			sigmai<- sigmai.n
			iter<- iter+1
		}
		u.mu<- 1/sqrt(sum(1/sigmai^2))
		N<- sum(n)
		llh<- -N*log(2*pi)/2-sum(log(sigmai^2))/2-sum((n*(x-mu)^2+(n-1)*si^2)/(2*sigmai^2))
		var.mu<- u.mu^2
	} else {
		# likelihood evaluated at the estimated parameters
		llh<- -N/2*log(2*pi)-sum(n*log(n))/2+sum(n*log(gammai/sigma2))/2-sum((n-1)*log(1-gammai))/2-sum(gammai*((x-mu)^2+(n-1)*s2/n/(1-gammai)))/sigma2/2
	}

	result<- list( mu=as.vector(mu), var.mu=as.vector(var.mu), sigma2=as.vector(sigma2), llh=as.vector(llh), tot.iter=as.vector(t), cur.rel.abs.error=as.vector(cur.rel.abs.error), gammai=as.vector(gammai), converged=converged, reduced.model=(sigma2==0) )

	class(result)<- "summary.vr.mle"

	return(result)

} # function block

print.summary.vr.mle <- function(x, ..., digits=3) {
	cat("\nVangel-Ruhkin estimate by MLE\n\n")

	est<-matrix(c(x$mu, x$var.mu, sqrt(x$var.mu)), ncol=3)
	dimnames(est) <- list("estimate", c("Value", "Var", "SD"))
	print(format(est, ..., digits=digits),quote=FALSE)

	cat(paste("\n Between-group variance: ", format(x$sigma2, ..., digits=digits)))
	cat(paste("\n         Log-likelihood: ", format(x$llh, ..., digits=digits)))
	cat(paste("\n             Iterations: ", format(x$tot.iter, ..., digits=digits)))
	cat(paste("\nRelative absolute error: ", format(x$cur.rel.abs.error, ..., digits=digits)))
	cat(paste("\n              Converged: ", format(x$converged, ..., digits=digits)))
	cat(paste("\n          Reduced model: ", format(x$reduced.model, ..., digits=digits)))
	
	
	cat("\n\nEstimated weights:\n")
	cat(paste(format(x$gammai, ..., digits=digits)))
	cat("\n\n")
}

