#TODO

.setUp <-function () {
	#library(MASS)
	.setUpDf <- within( list(),{
			x <- seq(0,1,length.out=41)[-c(1,41)]; #x[1] = x[1] + .Machine$double.eps; x[length(x)] <- x[length(x)]- .Machine$double.eps
			lx <- logit(x)			
		})
	attach(.setUpDf)
}

.tearDown <- function () {
	#detach(.setUpDf)
	detach()
}



test.inverseSame <- function(){
	xnorm <- logit(x)
	xinv <- invlogit(xnorm)
	checkEquals(x, xinv)
}

test.plogitnorm <- function(){
	px <- plogitnorm(x)	#percentiles
	checkEquals( pnorm(lx), px)
	px2 <- plogitnorm(x,mu=2,sigma=1)
	checkEquals( pnorm(lx,mean=2,sd=1), px2)
	#plot( px ~ x)
	#plot( px ~ logit(x))
	#lines( pnorm( logit(x)) ~ logit(x) )
}

test.twCoefLogitnorm <- function(){
	theta <- twCoefLogitnorm(0.7,0.9,perc=0.999)
	px <- plogitnorm(x,mu=theta[1],sigma=theta[2])	#percentiles function
	dx <- dlogitnorm(x,mu=theta[1],sigma=theta[2])	#density function
	#plot(px~x); abline(v=c(0.7,0.9)); abline(h=c(0.5,0.975))
	#plot(dx~x); abline(v=c(0.7,0.9))
	# upper percentile at 0.9
	checkEquals(which.min(abs(px-0.999)), which(x==0.9) )
	checkEquals(which.min(abs(px-0.5)), which.min(abs(x-0.7)) )
	# mode at 0.7
	#checkEquals(which(abs(x-0.7)<.Machine$double.eps), which.max(dx) )
}

test.twCoefLogitnormN <- function(){
	quant=c(0.7,0.8,0.9)
	perc=c(0.5,0.75,0.975)
	(theta <- twCoefLogitnormN( quant=quant, perc=perc ))
	#px <- plogitnorm(x,mu=theta[1],sigma=theta[2])	#percentiles function
	#dx <- dlogitnorm(x,mu=theta[1],sigma=theta[2])	#density function
	#plot(px~x); abline(v=quant,col="gray"); abline(h=perc,col="gray")
}


test.twCoefLogitnormMLE <- function(){
	theta <- twCoefLogitnormMLE(0.7,0.9,perc=0.975)
	px <- plogitnorm(x,mu=theta[1],sigma=theta[2])	#percentiles function
	dx <- dlogitnorm(x,mu=theta[1],sigma=theta[2])	#density function
	#plot(px~x); abline(v=c(0.7,0.9)); abline(h=c(0.5,0.975))
	#plot(dx~x); abline(v=c(0.7,0.9))
	# upper percentile at 0.9
	checkEquals(which.min(abs(px-0.975)), which(x==0.9) )
	# mode at 0.7
	checkEquals(which.min(abs(x-0.7)), which.max(dx) )
}


test.twCoefLogitnormE <- function(){
	theta <- twCoefLogitnormE(0.7,0.9)
	px <- plogitnorm(x,mu=theta[1],sigma=theta[2])	#percentiles function
	dx <- dlogitnorm(x,mu=theta[1],sigma=theta[2])	#density function
	#plot(px~x); abline(v=c(0.7,0.9)); abline(h=c(0.5,0.975))
	#plot(dx~x); abline(v=c(0.7,0.9))
	# upper percentile at 0.9
	checkEquals(which.min(abs(px-0.975)), which(x==0.9) )
	# mean at 0.7
	checkEqualsNumeric( momentsLogitnorm(mu=theta[1],sigma=theta[2])["mean"], 0.7, tolerance=1e-3)
	z <- rlogitnorm(1e5, mu=theta[1],sigma=theta[2])
	checkEqualsNumeric(0.7, mean(z), tolerance=5e-3 )
}

.tmp.f <- function(){
	px <- plogitnorm(x)	#percentiles
	plot( px ~ x )	
	plot( qlogitnorm(px) ~ x ) 	#one to one line
	plot( dlogitnorm(x,mu=0.9) ~ x, type="l" )
	abline( v=qlogitnorm(c(0.025,0.5,0.975), mu=0.9))
}

.tmp.f <- function(){
	library(MASS)
	?fitdistr	#not implemented
	quant = c(0.6,0.9)
	perc = c(0.5,0.975)
	theta0=c(mu=0,sigma=1)
	method="BFGS"
	#mtrace(ofLogitnorm)
	#popt <- as.list(tmp$par)
	popt <- as.list(coefLogitnorm(quant))
	popt2 <- as.list(coefLogitnorm(quant, perc=c(0.5,0.9995)))
	ofLogitnorm(popt,quant,perc)
	
	plot( dlogitnorm(x,mu=popt$mu,sigma=popt$sigma) ~ x, type="l" )
	abline( v=qlogitnorm(c(0.025,0.5,0.975),mu=popt$mu,sigma=popt$sigma))
	
	lines( dlogitnorm(x,mu=popt2$mu,sigma=popt2$sigma) ~ x, type="l", col="maroon" )
	abline( v=qlogitnorm(c(0.5,0.9995),mu=popt2$mu,sigma=popt2$sigma), col="maroon")
	
	popt <- as.list(coefLogitnorm(c(0.9, 0.9995), perc=c(0.5,0.9995)))
	plot( dlogitnorm(x,mu=popt$mu,sigma=popt$sigma) ~ x, type="l" )
	abline( v=qlogitnorm(c(0.025,0.5,0.975),mu=popt$mu,sigma=popt$sigma))
	
}

.tmp.f <- function(){
	#visualize the objective functions surface
	quant = c(0.6,0.9)
	perc = c(0.5,0.975)
	perc = c(0.5,0.9995)
	perc=c(0.5,upperBoundProb)
	quant = parms.var[varDist=="logitnorm",c("qMedian","qUpper")]
	quant = parms.var["epsF",c("qMedian","qUpper")]
	quant = parms.var["epsG",c("qMedian","qUpper")]
	quant = parms.var["epsP",c("qMedian","qUpper")]
	
	
	tmp.n <- 80
	tmp.mu <- seq(0,1,length.out=tmp.n)
	tmp.sigma <- seq(0.01,2.5,length.out=tmp.n)
	tmp <- as.matrix(expand.grid( mu=tmp.mu, sigma=tmp.sigma))
	tmp.of <- apply(tmp,1,ofLogitnorm, quant=quant,perc=perc )
	tmp.ofm <- matrix(tmp.of, nrow=tmp.n )
	image(tmp.mu, tmp.sigma, tmp.ofm)
	image(tmp.mu, tmp.sigma, exp(-0.5*tmp.ofm))	#very flat
	image(tmp.mu, tmp.sigma, exp(-0.5*1/(1/800)*tmp.ofm))  #distort by decreasing Temp *1/T
	tmp.o <- coefLogitnorm(quant=quant, perc=perc, returnDetails=TRUE)
	tmp.o
	popt <- as.list(tmp.o$par)
	points( popt$mu, popt$sigma)
	
	windows()
	unlist(popt)
	plot( dlogitnorm(x,mu=popt$mu,sigma=popt$sigma) ~ x, type="l" )
	abline( v=qlogitnorm(perc,mu=popt$mu,sigma=popt$sigma))
	lines( dlogitnorm(x,mu=popt$mu,sigma=1) ~ x, type="l", col="maroon" )
	abline( v=qlogitnorm(c(0.025,0.5,0.975),mu=popt$mu,sigma=1), col="maroon")
}

