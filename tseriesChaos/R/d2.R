#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:15:39 +0100 (ven, 02 dic 2005) $
d2 <- function(series, m, d, t, eps.min, neps=100) {
	checkEmbParms(series, m, d, t)
	if(eps.min<=0) stop("eps.min must be positive")
	neps <- as.integer(neps)
	if(neps<=0) neps <- 100
	res <- numeric(neps*m)
	eps.max <- diff(range(series))*sqrt(m)
	res <- .C("d2", series=as.double(series), length=as.integer(length(series)), m=as.integer(m), d=as.integer(d), t=as.integer(t), neps=as.integer(neps), eps.max=as.double(eps.max), eps.min=as.double(eps.min), res=as.double(res), 
	PACKAGE="tseriesChaos")[["res"]]
	res <- t(matrix(res, m, neps))
	res <- res[neps:1,]
	denom <- length(series) - (m-1)*d
	denom <- (denom-t+1)*(denom-t)/2
	res <- apply(res, 2, cumsum)/denom
	a <- -log(eps.min/eps.max)/(neps-1)
	eps <- eps.max*exp((1-1:neps)*a)
	eps <- eps[neps:1]
	res <- cbind(eps, res)
	colnames(res) <- c("eps",paste("m", 1:m,sep=""))
	class(res) <- "d2"
	res
}

plot.d2 <- function(x, ...) {
	m <- ncol(x)
	plot(x[,c(1,m)], type="l", log="xy", main="Sample correlation integral", xlab=expression(epsilon), ylab=expression(C(epsilon)))
	for(i in (m-1):2) lines(x[,c(1,i)])
}

print.d2 <- function(x, ...) {
	plot(x)
}
