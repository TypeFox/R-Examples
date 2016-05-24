normalmix.init=function (x, lambda = NULL, mu = NULL, s = NULL, k = 2, arbmean = TRUE, 
    arbvar = TRUE) 
{
    if (!is.null(s)) {
        arbvar <- (length(s) > 1)
        if (arbvar) 
            k <- length(s)
    }
    if (!is.null(mu)) {
        arbmean <- (length(mu) > 1)
        if (arbmean) {
            k <- length(mu)
            if (!is.null(s) && length(s) > 1 && k != length(s)) {
                stop("mu and sigma are each of length >1 but not of the same length.")
            }
        }
    }
    if (!arbmean && !arbvar) {
        stop("arbmean and arbvar cannot both be FALSE")
    }
    n = length(x)
    x = sort(x)
    x.bin = list()
    for (j in 1:k) {
        x.bin[[j]] <- x[max(1, floor((j - 1) * n/k)):ceiling(j * 
            n/k)]
    }
    if (is.null(s)) {
        s.hyp = as.vector(sapply(x.bin, sd))
	  if(any(s.hyp==0)) s.hyp[which(s.hyp==0)] = runif(sum(s.hyp==0),0,sd(x))
        if (arbvar) {
            s = 1/rexp(k, rate = s.hyp)
        }
        else {
            s = 1/rexp(1, rate = mean(s.hyp))
        }
    }
    if (is.null(mu)) {
        mu.hyp <- as.vector(sapply(x.bin, mean))
        if (arbmean) {
            mu = rnorm(k, mean = mu.hyp, sd = s)
        }
        else {
            mu = rnorm(1, mean = mean(mu.hyp), sd = mean(s))
        }
    }
    if (is.null(lambda)) {
        lambda <- runif(k)
        lambda <- lambda/sum(lambda)
    }
    else {
        lambda <- rep(lambda, length.out = k)
        lambda <- lambda/sum(lambda)
    }
    list(lambda = lambda, mu = mu, s = s, k = k, arbvar = arbvar, 
        arbmean = arbmean)
}


