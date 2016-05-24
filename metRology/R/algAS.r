#Algorithm A from ISO 5725-5:1998

algA <- function(x, k=1.5, na.rm=FALSE, tol=.Machine$double.eps^0.25,
	maxiter=25, verbose=FALSE) {

	if(na.rm) x <- na.omit(x)
	mu <- median(x)
	s <- ds <- mad(x)
	n <- length(x)
	if (s == 0) 
		stop("Initial estimated scale is zero.")
	
	theta <- 2*pnorm(abs(k))-1
	gamma <- 1/sqrt(theta+(1-theta)*k^2-2*k*dnorm(k))
		#gamma is identifiable as 1/sqrt(beta) where beta 
		#is the correction factor used on the variance by hubers

	iter <- 0
	if(verbose) cat(sprintf("%2d: mu=%f; s=%f\n", iter, mu, s))
	if(verbose>1) cat(paste(c("    z:", round(sort(x),2), "\n\n"), collapse=" "))
	while( (ds > tol*s) && iter<maxiter) {
		s.old <- s
		iter <- iter+1
		z <- pmin(pmax(mu - k * s, x), mu + k * s)
		mu <- mean(z)
		s <- gamma * sqrt(sum((z-mu)^2)/(n-1))
		ds <- abs(s.old-s)
		if(verbose) cat(sprintf("%2d: mu=%f; s=%f\n", iter, mu, s))
		if(verbose>1) cat(paste(c("    z:", round(sort(z),2), "\n\n"), collapse=" "))
	}

	if(iter >= maxiter) warning("Maximum iterations reached; algA may not have converged")

	return(list(mu=mu, s=s))
 	
}

#Algorithm S from ISO 5725-5:1998

algS <- function(s, degfree, na.rm=FALSE, prob.eta=0.9, is.range=FALSE, 
	tol=.Machine$double.eps^0.25, maxiter=25, verbose=FALSE) {
	#s is a vector of standard deviations
	#degfree is _the_ (scalar) number of degrees of freedom
	#associated with each
	#if degfree is a vector of length length(s), it will be 
	#replaced with median(degfree)
	
	#NA's in s are omitted if na.rm is TRUE
	
	#prob.eta is set to give the 10% upper tail of the chi-squared 
	#distribution for eta following Annex B of ISO 5725-5:1998. 
	
	#if is.range is TRUE, s is interpreted as a vector of ranges of 
	#two values and degfree is set to 1
	
	#tol is an arbitrary tolerance; iterations cease when the change
	#in estimated sd s.est is less than s.est*tol
	
	degfree <- if(!is.range) 
			median(degfree, na.rm=TRUE)
		   else
		   	1
		   	
	if(na.rm) s <- na.omit(s)
	
	p<-length(s)
	
	eta <- sqrt(qchisq(prob.eta, degfree)/degfree)
	
	z <- pchisq(degfree*eta*eta, degfree + 2, lower.tail=TRUE)
		
	xi <- 1/sqrt(z + (1-prob.eta)*eta*eta)
	
	s.est <- ds <- median(s)
	iter <- 0
	w<-s
	if(verbose) cat(sprintf("%2d: s.est=%f\n", iter, s.est))
	if(verbose>1) cat(paste(c("    w:", round(sort(w),2), "\n\n"), collapse=" "))
	while((ds > s.est*tol) && iter<maxiter) {
		iter <- iter+1
		s.old <- s.est
		psi <- eta * s.est
		w <- pmin(psi, s)
		s.est <- xi * sqrt( sum(w^2) / p )
		ds <- abs(s.old-s.est)
		if(verbose) cat(sprintf("%2d: psi=%f: s=%f\n", iter, psi, s.est))
		if(verbose>1) cat(paste(c("    w:", round(sort(w),2), "\n\n"), collapse=" "))
	}
	
	if(iter >= maxiter) warning("Maximum iterations reached; algS may not have converged")
	
	if(is.range) s.est <- s.est/sqrt(2)
	
	return(s.est)
}
