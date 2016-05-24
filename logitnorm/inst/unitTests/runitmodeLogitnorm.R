#TODO

.setUp <-function () {
}

.tearDown <- function () {
}

test.rightMode <- function(){
	theta0 <- c(mu=1.5, sigma=0.8)
	#plot the true and the rediscovered distributions
	xGrid = seq(0,1, length.out=81)[-c(1,81)]
	dx <- dlogitnorm(xGrid, mu=theta0[1], sigma=theta0[2])
	plot( dx~xGrid, type="l")
	
	mle <- modeLogitnorm(mu=theta0[1], sigma=theta0[2] )
	abline(v=mle,col="gray")
	
	#check by monte carlo integration
	#z <- rlogitnorm(1e6, mu=theta0[1], sigma=theta0[2]);	var(z)
	#dz <- density(z)
	#checkEqualsNumeric( dz$x[which.max(dz$y)], mle, tolerance=5e-2)
	checkEqualsNumeric( 0.88, mle, tolerance=1e-2)
}

test.leftMode <- function(){
	theta0 <- c(mu=-1.5, sigma=0.8)
	#plot the true and the rediscovered distributions
	xGrid = seq(0,1, length.out=81)[-c(1,81)]
	dx <- dlogitnorm(xGrid, mu=theta0[1], sigma=theta0[2])
	plot( dx~xGrid, type="l")
	
	mle <- modeLogitnorm(mu=theta0[1], sigma=theta0[2] )
	abline(v=mle,col="gray")
	
	#check by monte carlo integration
	# deprecated: did not run on Windows
	#z <- rlogitnorm(1e6, mu=theta0[1], sigma=theta0[2]);	var(z)
	#dz <- density(z)
	#checkEqualsNumeric( dz$x[which.max(dz$y)], mle, tolerance=5e-2)
	checkEqualsNumeric( 0.12, mle, tolerance=1e-2)		# regression 0.12 calculated previously
}

