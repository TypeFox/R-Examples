repnormmix.init = function (x, lambda = NULL, mu = NULL, s = NULL, k = 2, arbmean = TRUE, arbvar = TRUE) 
{
    n <- ncol(x)
    m <- nrow(x)
    y <- apply(x, 2, mean)
    x <- x[, order(y)]
    x.bin = list()
    for (j in 1:k) {
        x.bin[[j]] <- x[, max(1, floor((j - 1) * n/k)):ceiling(j * 
            n/k)]
    }
    if (is.null(s)) {
        s.hyp = sapply(lapply(x.bin, as.vector), sd)
        if (arbvar) {
            s = 1/rexp(k, rate = s.hyp)
        }
        else {
            s.hyp = mean(s.hyp)
            s = 1/rexp(1, rate = s.hyp)
        }
    }
    if (is.null(s) == FALSE && arbvar == TRUE) {
        k = length(s)
    }
    if (is.null(mu)) {
        mu.hyp <- sapply(lapply(x.bin, as.vector), mean)
        mu = rnorm(k, mean = mu.hyp, sd = s)
		if(arbmean==FALSE){
		mu = mean(mu)
		}
    }
	if (is.null(mu)==FALSE && arbmean==TRUE){
    		k = length(mu)
	}
    if (is.null(lambda)) {
        lambda = runif(k)
        lambda = lambda/sum(lambda)
    }
    else k = length(lambda)
    list(lambda = lambda, mu = mu, s = s, k = k)
}