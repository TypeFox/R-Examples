dpoislind <- function(x, theta, log = FALSE){
    if (theta <= 0) {
        stop(paste("theta must be positive!",
            "\n"))
    }
    p <- (theta^2*(x+theta+2)/(theta+1)^(x+3))*(x>=0)
    if (log)
        p <- log(p)
    p[is.nan(p)] <- 0
    if(!log) p <- pmin(pmax(p, 0),1)
    p
}

ppoislind <- function(q, theta, lower.tail = TRUE, log.p = FALSE){
    if (theta <= 0) {
        stop(paste("theta must be positive!",
            "\n"))
    }
	ind <- (q<0)
	q <- floor(q)
    temp <- sapply(1:length(q),function(i) sum(dpoislind(0:q[i],theta=theta,log=FALSE)))
    if (lower.tail == FALSE)
        temp <- 1 - temp
	if(any(ind)) temp[ind] <- 0 + 1*!lower.tail
    if (log.p)
        temp <- log(temp)
	if (!log.p)
		temp <- pmin(pmax(temp, 0),1)
	temp
}

qpoislind <- function(p, theta, lower.tail = TRUE, log.p = FALSE){
    if (theta <= 0) {
        stop(paste("theta must be positive!",
            "\n"))
    }
    if (log.p) p <- exp(p)

	up <- ifelse(theta>.125,400,2000)
	if(lower.tail){
		tmp <- ppoislind(0:up,theta=theta)
		all.p <- suppressWarnings(sapply(1:length(p),function(i) min(which(tmp>=p[i]))-1))
		all.p[p==1] <- Inf
		all.p[p==0] <- 0
		all.p[(p>1)|(p<0)] <- NaN
	} else{
		tmp <- ppoislind(0:up,theta=theta,lower.tail=FALSE)
		all.p <- suppressWarnings(pmax(sapply(1:length(p),function(i) max(which(tmp>p[i]))),0))
		if(up==2000&any(all.p==2000)) all.p[all.p==2000] <- Inf 
		all.p[p==1] <- 0
		all.p[p==0] <- Inf
		all.p[(p>1)|(p<0)] <- NaN
	}
	if(any(is.nan(all.p))) warning("NaNs produced")
	all.p
}

rpoislind <- function(n, theta){
    if (theta <= 0) {
        stop(paste("theta must be positive!",
            "\n"))
    }
	u <- runif(n)
	p <- theta/(theta+1)
	ind <- (u>p)
	lambda <- rexp(n,theta)+(rexp(n,theta))*ind
	out <- rpois(n,lambda)
	out
}

