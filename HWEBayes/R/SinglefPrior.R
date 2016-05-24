SinglefPrior <-
function(nsim,alpha,lambdamu,lambdasd){
	    p <- rdirichlet(n=nsim,alpha=alpha)
	    lambda <- rnorm(nsim,mean=lambdamu,sd=lambdasd)
	    lgts <- apply(p, 1, function(p) {baselogit(p)$baselogit})
            if (is.matrix(lgts)) lgts <- t(lgts)
	    pmin <- apply(p,1,min)
	    fmin <- -pmin/(1-pmin)
	    f <- (exp(lambda)+fmin)/(exp(lambda)+1)
	    list(p=p,f=f,lgts=lgts,lambda=lambda)
}

