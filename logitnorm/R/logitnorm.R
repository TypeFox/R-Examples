# random, percentile, density and quantile function of the logit-normal distribution
# and estimation of parameters from percentiles by Sum of Squares Newton optimization
# 

logit <- function(
	### Transforming (0,1) to normal scale (-Inf Inf)
	p,...
){ 
	##details<<
	## function \eqn{ logit(p)= log \left( \frac{p}{1-p} \right) = log(p) - log(1-p) }
	##seealso<< \code{\link{invlogit}}
	##seealso<< \code{\link{logitnorm}}
	qlogis(p,...) 
}

invlogit <- function(
	### Transforming (-Inf,Inf) to original scale (0,1)
	q,...
){
	##details<<
	## function \eqn{f(z) = \frac{e^{z}}{e^{z} + 1} \! = \frac{1}{1 + e^{-z}} \!}
	##seealso<< \code{\link{logit}}
	##seealso<< \code{\link{logitnorm}}
	plogis(q,...) 
} 

# for normal-logit transformation do not use scale and location, better use it in generating rnorm etc

rlogitnorm <- function(
	### Random number generation for logitnormal distribution
	mu = 0, sigma = 1,	##<< distribution parameters
	...	##<< arguments to \code{\link{rnorm}}
){	
	##seealso<< \code{\link{logitnorm}}
	plogis( rnorm(mean=mu,sd=sigma,...) ) 
}

plogitnorm <- function(
	### Distribution function for logitnormal distribution	
	q
	,mu = 0, sigma = 1	##<< distribution parameters
	,...
){
	##seealso<< \code{\link{logitnorm}}
	ql <- qlogis(q)
	pnorm(ql,mean=mu,sd=sigma,...)
}

dlogitnorm <- function(
	### Density function of logitnormal distribution	
	q		##<< quantiles
	,mu = 0, sigma = 1	##<< distribution parameters
	,...	##<< further arguments passed to \code{\link{dnorm}}: \code{mean}, and \code{sd} for mu and sigma respectively.  
){
	##alias<< logitnorm
	
	##details<< \describe{\item{Logitnorm distribution}{ 
	## \itemize{
	## \item{density function: dlogitnorm }
	## \item{distribution function: \code{\link{plogitnorm}} }
	## \item{quantile function: \code{\link{qlogitnorm}} }
	## \item{random generation function: \code{\link{rlogitnorm}} }
	## }
	## }}
	
	##details<< \describe{\item{Transformation functions}{ 
	## \itemize{
	## \item{ (0,1) -> (-Inf,Inf): \code{\link{logit}} }
	## \item{ (-Inf,Inf) -> (0,1): \code{\link{invlogit}} }
	## }
	## }}
	
	##details<< \describe{\item{Moments and mode}{ 
	## \itemize{
	## \item{ Expected value and variance: \code{\link{momentsLogitnorm}} }
	## \item{ Mode: \code{\link{modeLogitnorm}} }
	## }
	## }}
	
	##details<< \describe{\item{Estimating parameters}{ 
	## \itemize{
	## \item{from mode and upper quantile: \code{\link{twCoefLogitnormMLE}} }
	## \item{from median and upper quantile: \code{\link{twCoefLogitnorm}} }
	## \item{from expected value, i.e. mean and upper quantile: \code{\link{twCoefLogitnormE}} }
	## \item{from a confidence interval which is symmetric at normal scale: \code{\link{twCoefLogitnormCi}} }
	## \item{from prescribed quantiles: \code{\link{twCoefLogitnormN}} }
	## }
	## }}
	
	ql <- qlogis(q)
	dnorm(ql,mean=mu,sd=sigma,...) /q/(1-q)	# multiply by the Jacobian (derivative) of the back-Transformation (logit)
}

qlogitnorm <- function(
	### Quantiles of logitnormal distribution.	
	p
	,mu = 0, sigma = 1	##<< distribution parameters
	,...
){
	##seealso<< \code{\link{logitnorm}}
	qn <- qnorm(p,mean=mu,sd=sigma,...)
	plogis(qn)
}


#------------------ estimating parameters from fitting to distribution 
.ofLogitnorm <- function(
	### Objective function used by \code{\link{coefLogitnorm}}. 
	theta				##<< theta[1] is mu, theta[2] is the sigma
	,quant				##<< q are the quantiles for perc
	,perc=c(0.5,0.975) 
){
	# there is an analytical expression for the gradient, but here just use numerical gradient
	
	# calculating perc should be faster than calculating quantiles of the normal distr.
	if( theta[2] <= 0) return( Inf )
	tmp.predp = plogitnorm(quant, mu=theta[1], sigma=theta[2] )
	tmp.diff = tmp.predp - perc
	
	#however, q not so long valley and less calls to objective function
	#but q has problems with high sigma
	#tmp.predq = qlogitnorm(perc, mean=theta[1], sigma=theta[2] )
	#tmp.diff =  tmp.predq - quant
	sum(tmp.diff^2)
}

twCoefLogitnormN <- function(
	### Estimating coefficients of logitnormal distribution from a vector of quantiles and perentiles (non-vectorized).	
	quant					##<< the quantile values
	,perc=c(0.5,0.975)		##<< the probabilites for which the quantiles were specified
	,method="BFGS"			##<< method of optimization (see \code{\link{optim}})
	,theta0=c(mu=0,sigma=1)	##<< starting parameters
	,returnDetails=FALSE	##<< if TRUE, the full output of optim is returned instead of only entry par
	, ... 					##<< further parameters passed to optim, e.g. control=list(maxit=1000)
){
	##seealso<< \code{\link{logitnorm}}
	names(theta0) <- c("mu","sigma")
	tmp <- optim( theta0, .ofLogitnorm, quant=quant,perc=perc, method=method, ...)
	if( tmp$convergence != 0)
		warning(paste("coefLogitnorm: optim did not converge. theta0=",paste(theta0,collapse=",")))
	if( returnDetails )
		tmp
	else
		tmp$par
	### named numeric vector with estimated parameters of the logitnormal distrubtion.
	### names: \code{c("mu","sigma")}
}
#mtrace(coefLogitnorm)
attr(twCoefLogitnormN,"ex") <- function(){
	# estimate the parameters
	quant=c(0.7,0.8,0.9)
	perc=c(0.5,0.75,0.975)
	(theta <- twCoefLogitnormN( quant=quant, perc=perc ))
	
	x <- seq(0,1,length.out=41)[-c(1,41)]	# plotting grid
	px <- plogitnorm(x,mu=theta[1],sigma=theta[2])	#percentiles function
	plot(px~x); abline(v=quant,col="gray"); abline(h=perc,col="gray")
}

twCoefLogitnorm <- function(
	### Estimating coefficients of logitnormal distribution from median and upper quantile	
	median					##<< numeric vector: the median of the density function
	,quant					##<< numeric vector: the upper quantile value
	,perc=0.975				##<< numeric vector: the probability for which the quantile was specified
	,method="BFGS"			##<< method of optimization (see \code{\link{optim}})
	,theta0=c(mu=0,sigma=1)	##<< starting parameters
	,returnDetails=FALSE	##<< if TRUE, the full output of optim is attached as attributes resOptim
	, ... 
){
	##seealso<< \code{\link{logitnorm}}
	# twCoefLogitnorm
	names(theta0) <- c("mu","sigma")
	nc <- c(length(median),length(quant),length(perc)) 
	n <- max(nc)
	res <- matrix( as.numeric(NA), n,2, dimnames=list(NULL,c("mu","sigma")))
	resOptim <- list()
	qmat <- cbind(median,quant)
	pmat <- cbind(0.5,perc)
	for( i in 1:n ){	# for each row in (recycled) vector
		i0 <- i-1
		tmp <- optim( theta0, .ofLogitnorm, perc=pmatI<-as.numeric(pmat[1+i0%%nrow(pmat),]), quant=(qmatI<-qmat[1+i0%%nrow(qmat),]), method=method, ...)
		if( tmp$convergence == 0)
			res[i,] <- tmp$par
		else
			warning(paste("coefLogitnorm: optim did not converge. theta0=",paste(theta0,collapse=","),"; median=",qmatI[1],"; quant=",qmatI[2],"; perc=",pmatI[2],sep=""))
		resOptim[[i]] <- tmp
	}
	if( returnDetails ) attr(res,"resOptim") <- resOptim
	res
	### numeric matrix with columns \code{c("mu","sigma")}
	### rows correspond to rows in median, quant, and perc
}
#mtrace(twCoefLogitnorm)
attr(twCoefLogitnorm,"ex") <- function(){
	# estimate the parameters, with median at 0.7 and upper quantile at 0.9
	(theta <- twCoefLogitnorm(0.7,0.9))
	
	x <- seq(0,1,length.out=41)[-c(1,41)]	# plotting grid
	px <- plogitnorm(x,mu=theta[1],sigma=theta[2])	#percentiles function
	plot(px~x); abline(v=c(0.7,0.9),col="gray"); abline(h=c(0.5,0.975),col="gray")
	
	dx <- dlogitnorm(x,mu=theta[1],sigma=theta[2])	#density function
	plot(dx~x); abline(v=c(0.7,0.9),col="gray")
	
	# vectorized
	(theta <- twCoefLogitnorm(seq(0.4,0.8,by=0.1),0.9))
}

twCoefLogitnormCi <- function( 
	### Calculates mu and sigma of the logitnormal distribution from lower and upper quantile, i.e. confidence interval.
	lower	##<< value at the lower quantile, i.e. practical minimum
	,upper	##<< value at the upper quantile, i.e. practical maximum
	,perc=0.975	##<< numeric vector: the probability for which the quantile was specified
	,sigmaFac=qnorm(perc) 	##<< sigmaFac=2 is 95% sigmaFac=2.6 is 99% interval
	,isTransScale = FALSE ##<< if true lower and upper are already on logit scale
){	
	##seealso<< \code{\link{logitnorm}}
	if( !isTRUE(isTransScale) ){
		lower <- logit(lower)
		upper <- logit(upper)
	}
	halfWidth <- (upper-lower)/2
	sigma <- halfWidth/sigmaFac
	cbind( mu=upper-halfWidth, sigma=sigma )
	### named numeric vector: mu and sigma parameter of the logitnormal distribution.
}
attr(twCoefLogitnormCi,"ex") <- function(){
	mu=2
	sd=c(1,0.8)
	p=0.99
	lower <- l <- qlogitnorm(1-p, mu, sd )		# p-confidence interval
	upper <- u <- qlogitnorm(p, mu, sd )		# p-confidence interval
	cf <- twCoefLogitnormCi(lower,upper)	
	all.equal( cf[,"mu"] , c(mu,mu) )
	all.equal( cf[,"sigma"] , sd )
}

.ofLogitnormMLE <- function(
	### Objective function used by \code{\link{coefLogitnormMLE}}. 
	mu	##<< numeric vector of proposed parameter				
		## make sure that logit(mu)<mle for mle>0.5 and logit(mu)>mle for mle<0.5
	,mle				##<< the mode of the density distribution
	,logitMle=logit(mle)##<< may provide for performance reasons
	,quant				##<< q are the quantiles for perc
	,perc=c(0.975) 
){
	# given mu and mle, we can calculate sigma
	sigma2 = (logitMle-mu)/(2*mle-1)
	ifelse( sigma2 <= 0, Inf, {
		tmp.predp = plogitnorm(quant, mu=mu, sigma=sqrt(sigma2) )
		tmp.diff = tmp.predp - perc
		tmp.diff^2
	})
}

twCoefLogitnormMLE <- function(
	### Estimating coefficients of logitnormal distribution from mode and upper quantile	
	mle						##<< numeric vector: the mode of the density function
	,quant					##<< numeric vector: the upper quantile value
	,perc=0.999				##<< numeric vector: the probability for which the quantile was specified
	#, ... 					##<< Further parameters to \code{\link{optimize}}
){
	##seealso<< \code{\link{logitnorm}}
	# twCoefLogitnormMLE
	nc <- c(length(mle),length(quant),length(perc)) 
	n <- max(nc)
	res <- matrix( as.numeric(NA), n,2, dimnames=list(NULL,c("mu","sigma")))
	for( i in 1:n ){
		i0 <- i-1
		mleI<-mle[1+i0%%nc[1]]
		if( mleI==0.5) mleI=0.5-.Machine$double.eps	# in order to estimate sigma
		# we now that mu is in (logit(mle),0) for mle < 0.5 and in (0,logit(mle)) for mle > 0.5 for unimodal distribution
		# there might be a maximum in the middle and optimized misses the low part
		# hence, first get near the global minimum by a evaluating the cost at a grid
		# grid is spaced narrower at the edge
		logitMle <- logit(mleI)
		#intv <- if( logitMle < 0) c(logitMle+.Machine$double.eps,0) else c(0,logitMle-.Machine$double.eps)
		upper <- abs(logitMle)-.Machine$double.eps
		tmp <- sign(logitMle)*log(seq(1,exp(upper),length.out=40))
		oftmp<-.ofLogitnormMLE(tmp, mle=(mleI), logitMle=logitMle, quant=(quantI<-quant[1+i0%%nc[2]]), perc=(percI<-perc[1+i0%%nc[3]]))
		#plot(tmp,oftmp)
		imin <- which.min(oftmp)
		intv <- tmp[ c(max(1,imin-1), min(length(tmp),imin+1)) ]
		if( diff(intv)==0 ) mu <- intv[1] else
			mu <- optimize( .ofLogitnormMLE, interval=intv, mle=(mleI), logitMle=logitMle, quant=(quantI<-quant[1+i0%%nc[2]]), perc=(percI<-perc[1+i0%%nc[3]]))$minimum
		sigma <- sqrt((logitMle-mu)/(2*mleI-1))
		res[i,] <- c(mu,sigma)
	}
	res
	### numeric matrix with columns \code{c("mu","sigma")}
	### rows correspond to rows in mle, quant, and perc
}
#mtrace(coefLogitnorm)
attr(twCoefLogitnormMLE,"ex") <- function(){
	
	# estimate the parameters, with mode 0.7 and upper quantile 0.9
	(theta <- twCoefLogitnormMLE(0.7,0.9))
	
	x <- seq(0,1,length.out=41)[-c(1,41)]	# plotting grid
	px <- plogitnorm(x,mu=theta[1],sigma=theta[2])	#percentiles function
	plot(px~x); abline(v=c(0.7,0.9),col="gray"); abline(h=c(0.999),col="gray")
	dx <- dlogitnorm(x,mu=theta[1],sigma=theta[2])	#density function
	plot(dx~x); abline(v=c(0.7,0.9),col="gray")
	
	# vectorized
	(theta <- twCoefLogitnormMLE(mle=seq(0.4,0.8,by=0.1),quant=0.9))
}

.ofLogitnormE <- function(
	### Objective function used by \code{\link{coefLogitnormE}}. 
	theta				##<< theta[1] is the mu, theta[2] is the sigma
	,mean				##<< the expected value of the density distribution
	,quant				##<< q are the quantiles for perc
	,perc=c(0.975) 
	,...				##<< further argument to \code{\link{momentsLogitnorm}}.
){
	tmp.predp = plogitnorm(quant, mu=theta[1], sigma=theta[2] )
	tmp.diff = tmp.predp - perc
	.exp <- momentsLogitnorm(theta[1],theta[2],...)["mean"]
	tmp.diff.e <- mean-.exp
	sum(tmp.diff^2) + tmp.diff.e^2
}


twCoefLogitnormE <- function(
	### Estimating coefficients of logitnormal distribution from expected value, i.e. mean, and upper quantile.	
	mean					##<< the expected value of the density function
	,quant					##<< the quantile values
	,perc=c(0.975)			##<< the probabilites for which the quantiles were specified
	,method="BFGS"			##<< method of optimization (see \code{\link{optim}})
	,theta0=c(mu=0,sigma=1)	##<< starting parameters
	,returnDetails=FALSE	##<< if TRUE, the full output of optim is returned with attribut resOptim
	, ... 
){
	##seealso<< \code{\link{logitnorm}}
	# twCoefLogitnormE
	names(theta0) <- c("mu","sigma")
	mqp <- cbind(mean,quant,perc)
	n <- nrow(mqp)
	res <- matrix( as.numeric(NA), n,2, dimnames=list(NULL,c("mu","sigma")))
	resOptim <- list()
	for( i in 1:n ){
		mqpi<- mqp[i,]
		tmp <- optim( theta0, .ofLogitnormE, mean=mqpi[1], quant=mqpi[2], perc=mqpi[3], method=method, ...)
		resOptim[[i]] <- tmp
		if( tmp$convergence != 0)
			warning(paste("coefLogitnorm: optim did not converge. theta0=",paste(theta0,collapse=",")))
		else
			res[i,] <- tmp$par
	}
	if( returnDetails )
		attr(res,"resOptim") <- resOptim
	res
	### named numeric matrix with estimated parameters of the logitnormal distrubtion.
	### colnames: \code{c("mu","sigma")}
}
#mtrace(coefLogitnorm)
attr(twCoefLogitnormE,"ex") <- function(){
	# estimate the parameters
	(thetaE <- twCoefLogitnormE(0.7,0.9))
	
	x <- seq(0,1,length.out=41)[-c(1,41)]	# plotting grid
	px <- plogitnorm(x,mu=thetaE[1],sigma=thetaE[2])	#percentiles function
	plot(px~x); abline(v=c(0.7,0.9),col="gray"); abline(h=c(0.5,0.975),col="gray")
	dx <- dlogitnorm(x,mu=thetaE[1],sigma=thetaE[2])	#density function
	plot(dx~x); abline(v=c(0.7,0.9),col="gray")
	
	z <- rlogitnorm(1e5, mu=thetaE[1],sigma=thetaE[2])
	mean(z)	# about 0.7
	
	# vectorized
	(theta <- twCoefLogitnormE(mean=seq(0.4,0.8,by=0.1),quant=0.9))
}

momentsLogitnorm <- function(
	### First two moments of the logitnormal distribution by numerical integration
	mu	##<< parameter mu 
	,sigma		##<< parameter sigma
	,abs.tol=0	##<< chaning default to \code{\link{integrate}}
	,...		##<< further parameters to the \code{\link{integrate}} function
){
	fExp <- function(x)  plogis(x)*dnorm(x,mean=mu,sd=sigma)
	.exp <- integrate(fExp,-Inf,Inf, abs.tol=abs.tol, ...)$value
	fVar <- function(x)   (plogis(x) - .exp)^2 * dnorm(x,mean=mu,sd=sigma)
	.var <- integrate(fVar,-Inf,Inf, abs.tol=abs.tol, ...)$value
	c( mean=.exp, var=.var )
	### named numeric vector with components \itemize{
	### \item{ \code{mean}: expected value, i.e. first moment}
	### \item{ \code{var}: variance, i.e. second moment }
	### }
}
attr(momentsLogitnorm,"ex") <- function(){
	(res <- momentsLogitnorm(4,1))
	(res <- momentsLogitnorm(5,0.1))
}

.ofModeLogitnorm <- function(
	### Objective function used by \code{\link{coefLogitnormMLE}}.
	mle		##<< proposed mode
	,mu	##<< parameter mu 
	,sd2	##<< parameter sigma^2
){
	if( mle <=0 | mle >=1) return(.Machine$double.xmax)
	tmp.diff.mle <- mu-sd2*(1-2*mle)   - logit(mle)
	tmp.diff.mle^2
}

modeLogitnorm <- function(
	### Mode of the logitnormal distribution by numerical optimization
	mu	##<< parameter mu 
	,sigma		##<< parameter sigma
	,tol=invlogit(mu)/1000	##<< precisions of the estimate
){
	##seealso<< \code{\link{logitnorm}}
	itv <- if(mu<0) c(0,invlogit(mu)) else c(invlogit(mu),1)
	tmp <- optimize( .ofModeLogitnorm, interval=itv, mu=mu, sd2=sigma^2, tol=tol)
	tmp$minimum
}

twSigmaLogitnorm <- function(
	### Estimating coefficients of logitnormal distribution from mode and given mu	
	mle		##<< numeric vector: the mode of the density function
	,mu=0	##<< for mu=0 the distribution will be the flattest case (maybe bimodal) 
){
	sigma = sqrt( (logit(mle)-mu)/(2*mle-1) )
	##seealso<< \code{\link{logitnorm}}
	# twSigmaLogitnorm
	cbind(mu=mu, sigma=sigma)
	### numeric matrix with columns \code{c("mu","sigma")}
	### rows correspond to rows in mle and mu
}
#mtrace(coefLogitnorm)
attr(twSigmaLogitnorm,"ex") <- function(){
	
	# estimate the parameters, with mode 0.7 and upper quantile 0.9
	(theta <- twCoefLogitnormMLE(0.7,0.9))
	
	x <- seq(0,1,length.out=41)[-c(1,41)]	# plotting grid
	px <- plogitnorm(x,mu=theta[1],sigma=theta[2])	#percentiles function
	plot(px~x); abline(v=c(0.7,0.9),col="gray"); abline(h=c(0.999),col="gray")
	dx <- dlogitnorm(x,mu=theta[1],sigma=theta[2])	#density function
	plot(dx~x); abline(v=c(0.7,0.9),col="gray")
	
	# vectorized
	(theta <- twCoefLogitnormMLE(mle=seq(0.4,0.8,by=0.1),quant=0.9))
}





