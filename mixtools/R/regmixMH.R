regmixMH=function (y, x, lambda = NULL, beta = NULL, s = NULL, k = 2, 
    addintercept = TRUE, mu = NULL, sig = NULL, lam.hyp = NULL, sampsize = 1000, 
    omega = 0.01, thin = 1) 
{
    if (addintercept) {
        x = cbind(1, x)
    }
    XTX <- solve(t(x)%*%x) 
    n <- length(y)
    p <- ncol(x)
    if (is.null(s)) {
        s = sqrt(1/rexp(k))
    }
    else k = length(s)
    if (is.null(beta)) {
        beta = matrix(rnorm(p * k), p, k)
    }
    else k = ncol(beta)
    if (is.null(lambda)) {
        lambda = runif(k)
        lambda = lambda/sum(lambda)
    }
    else k = length(lambda)
    if (is.null(mu)) {
        mu = 0 * beta
    }
    if (is.null(sig)) {
        sig = rep(5 * sqrt(var(y)),k)
    }
    sig.beta = t(matrix(rep(sig,p),ncol=p)) * sqrt(matrix(rep(diag(XTX),k),ncol=k)) 
    if(is.null(lam.hyp)) lam.hyp = rep(1,k)
    L.theta <- matrix(nrow = n, ncol = k)
    pi.beta <- matrix(nrow = p, ncol = k)
    pi.sigma <- c()
    pi.lambda <- c()
    new.L.theta <- matrix(nrow = length(y), ncol = k)
    new.pi.beta <- matrix(nrow = p, ncol = k)
    new.pi.sigma <- c()
    new.pi.lambda <- c()
    accepts = 0
    theta <- matrix(c(beta, s, lambda), nrow = 1)
    thetalist <- matrix(theta, sampsize, ncol(theta), byrow=TRUE)
    for (i in 2:sampsize) {
        log.pi.beta <- dnorm(beta, mu, sig.beta, log=TRUE)
        log.pi.sigma <- dexp(s, 1/sig, log=TRUE)
        log.pi.lambda <- sum((lam.hyp-1)*log(lambda)) - lgamma(lam.hyp) +
                         lgamma(sum(lam.hyp)) # Dirichlet log-density
#        log.pi.lambda <- log(ddirichlet(lambda, lam.hyp)) 
        L.theta <- dnorm(y - x %*% beta, 0, matrix(s, n, k, byrow = TRUE)
                         ) %*% matrix(lambda, k, 1)
        log.Lik.theta <- sum(log(L.theta))
        log.prior <- sum(log.pi.beta) + sum(log.pi.sigma) + sum(log.pi.lambda)
        log.f.theta <- log.Lik.theta + log.prior
        new.beta <- beta + omega * matrix(rcauchy(k * p), p, k)
        new.sigma <- log(s) + omega * rcauchy(k) 
        new.sigma <- exp(new.sigma) 
	  new.lambda <- lambda.pert(lambda, omega*rcauchy(k)) 
        log.new.pi.beta <- dnorm(new.beta, mu, sig.beta, log=TRUE) 
        log.new.pi.sigma <- dexp(new.sigma, 1/sig, log=TRUE) 
        log.new.pi.lambda <- sum((lam.hyp-1)*log(new.lambda)) - lgamma(lam.hyp) +
                              lgamma(sum(lam.hyp)) # Dirichlet log-density
#        log.new.pi.lambda <- log(ddirichlet(new.lambda, lam.hyp)) 
        new.L.theta <- dnorm(y - x %*% new.beta, 0, matrix(new.sigma, 
            n, k, byrow = TRUE)) %*% matrix(new.lambda, k, 1)
        log.new.Lik.theta <- sum(log(new.L.theta))
        log.new.prior <- sum(log.new.pi.beta) + sum(log.new.pi.sigma) + 
            sum(log.new.pi.lambda)
        log.new.f.theta <- log.new.Lik.theta + log.new.prior
        new.theta <- c(new.beta, new.sigma, new.lambda)
        a <- log.new.f.theta - log.f.theta
        r <- log(runif(1))
        if (a > r & !is.na(a)) {
            theta <- new.theta
            beta <- new.beta
            s <- new.sigma
            lambda <- new.lambda
            accepts = accepts + 1
        }
        if (i%%thin == 0) 
            thetalist[i,] = theta
    }
    cat(paste("Acceptance rate: ", 100 * accepts/sampsize, "%\n", 
        sep = ""))
    colnames(thetalist) <- c(paste("beta", rep(0:(p - 1), k), 
        ".", rep(1:k, rep(p, k)), sep = ""), paste("s.", 1:k, 
        sep = ""), paste("lambda.", 1:k, sep = ""))
    invisible(thetalist)
	a=list(x=x, y=y, theta=thetalist, components=k)
	class(a)="mixMCMC"
	a
}
