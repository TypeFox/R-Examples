"scenNID" <-
function(y, n, cz, type="L")
{
	m <- length(y)
	nc <- n-m
	ip <- x1 <- x2 <- numeric(n)
	if (nc > 0) {
	}
	ip[m+1:n] <- 1
	x1[1:m] <- y
	x2[(m+1):n] <- cz
	xmean <- mean(y)
	xsigma <- sd(y)
	e1 <- e2 <- 10^(-6)
	cov <- matrix(numeric(4), ncol=2, nrow=2)
	maxits <- 100
	k <- 0
  nobs <- c(0,0,0,0)
	ifault <- 0
	outF <- .Fortran("em", as.integer(n), as.double(x1), as.double(x2), 
			as.integer(ip), as.double(xmean), as.double(xsigma), 
			as.double(e1), as.double(e2), as.integer(maxits),
			as.double(cov), as.integer(nobs), as.integer(k), as.integer(ifault),
			PACKAGE="cents")
    muHat <- outF[[5]]
	sigHat <- outF[[6]]
	sds <- sqrt((outF[[10]])[c(1,4)])
	est <- matrix(c(muHat,sigHat,sds), ncol=2, nrow=2)
	dimnames(est) <- list(c("mean", "sd"), c("mle", "se(mle)"))
  ans <- list(est=est, nobs=outF[[11]], iterCount=outF[[12]], ifault=outF[[13]])
  ans
}

