mvnormalmix.init = function (x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, arbmean = TRUE, arbvar = TRUE) 
{
    n <- nrow(x)
    p <- ncol(x)
    y <- apply(x, 1, mean)
    x <- x[order(y), ]
    x.bin <- list()
    for (j in 1:k) {
        x.bin[[j]] <- x[max(1, floor((j - 1) * n/k)):ceiling(j * 
            n/k), ]
    }
    if (is.null(sigma)) {
        if (arbvar) {
            sigma.hyp = lapply(1:k, function(i) (apply(x.bin[[i]], 
                2, var))^-1)
            sigma = lapply(1:k, function(i) diag(1/rexp(p, rate = sigma.hyp[[i]])))
        }
        else {
            sigma.hyp = apply(sapply(1:k, function(i) (apply(x.bin[[i]], 
                2, var))^-1), 2, mean)
            sigma = diag(1/rexp(p, rate = sigma.hyp))
        }
    }
    if (is.null(sigma) == FALSE && arbvar == TRUE) {
        k = length(sigma)
    }
    if (is.null(mu)) {
        mu.hyp <- lapply(1:k, function(i) apply(x.bin[[i]], 2, 
            mean))
        if (arbvar) {
            mu <- lapply(1:k, function(i) as.vector(rmvnorm(1, 
                mu = as.vector(mu.hyp[[i]]), sigma = as.matrix(sigma[[i]]))))
        }
        else mu <- lapply(1:k, function(i) as.vector(rmvnorm(1, 
            mu = as.vector(mu.hyp[[i]]), sigma = as.matrix(sigma))))
	  if (arbmean==FALSE) {
		mu <- apply(sapply(mu,as.vector),1,mean)
#		mu <- lapply(1:k, function(i) mu)
    }
	}
    if (is.null(mu) == FALSE && arbmean == TRUE){
	k=length(mu)
	}
    if (is.null(lambda)) {
        lambda <- runif(k)
        lambda <- lambda/sum(lambda)
    }
    else k <- length(lambda)
    list(lambda = lambda, mu = mu, sigma = sigma, k = k)
}