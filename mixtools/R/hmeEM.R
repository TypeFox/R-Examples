hmeEM <- function (y, x, lambda=NULL, beta = NULL, sigma = NULL, w = NULL, 
        k = 2, addintercept = TRUE, epsilon = 1e-08, maxit = 10000, verb=FALSE)
    {
    NR <- function (w, x, z)
    {
    A = (inv.logit(x %*% w)/(1 + exp(x %*% w)))^-1
    B = z[,1]-inv.logit(x %*% w)
    C = matrix(c(sum(B),sum(x[,2]*B)),2,1)
    D = matrix(nrow=2,ncol=nrow(x))
    D = solve(matrix(c(sum(A),sum(A * x[,2]),sum(A * x[,2]),sum(A * x[,2]^2)),nrow=2,ncol=2))
    E = w + apply(D %*% C, 1, sum)
    list(w.new = E)
    }
    if (k != 2) cat("This algorithm currently only works for a 2-component mixture-of-experts model!")
    s=sigma
    if (addintercept) {
        x = cbind(1, x)
    }
    n <- length(y)
    p <- ncol(x)
    if (is.null(w)) {
        w = rep(0,p*(k-1))
    }

    tmp <- regmix.init(y = y, x = x, beta = beta, 
        s = s, k = k, addintercept = addintercept)
    beta = tmp$beta
    s = tmp$s
    k = tmp$k    
    if(is.null(lambda)) lambda = inv.logit(x %*% w)
    diff <- 1
    iter <- 0
    nlogroot2pi = n * log(sqrt(2 * pi))
    xbeta <- x %*% beta
    tmp <- t(-1/2/s[2]^2 * t(y - xbeta[, -1])^2) + 1/2/s[1]^2 *
        (y - xbeta[, 1])^2
    tmp2 <- (1-lambda)/s[2] * s[1]/lambda
    tmp3 <- log(1 + (tmp2 * exp(tmp)))
    obsloglik <- sum(log(lambda)) - nlogroot2pi - n * log(s[1]) -
        1/2/s[1]^2 * sum((y - xbeta[, 1])^2) + sum(tmp3)
    ll <- obsloglik
    z = matrix(nrow = n, ncol = k)
	restarts <- 0
    while (diff > epsilon & iter < maxit) {

    lambda.mat <- cbind(lambda,1-lambda)

    res <- (y - xbeta)^2

        for (i in 1:n) {
            for (j in 1:k) {
                z.denom = c()
                for (h in 1:k) {
                  z.denom = c(z.denom, (lambda.mat[i,h]/lambda.mat[i,j]) * 
                    (s[j]/s[h]) * exp(-0.5 * ((1/s[h]^2) * res[i, h] - (1/s[j]^2) * res[i, j])))
                }
                z[i, j] = 1/sum(z.denom)
            }
        }
#        z[, k] = 1 - apply(as.matrix(z[, (1:(k - 1))]), 1, sum)
	z = z/apply(z,1,sum)

                if (addintercept) {
                  lm.out <- lapply(1:k, function(i) lm(y ~ x[, 
                    -1], weights = z[, i]))
                }
                else lm.out <- lapply(1:k, function(i) lm(y ~ 
                  x - 1, weights = z[, i]))
                beta.new <- sapply(lm.out, coef)


	w.diff=10
	while(sum(abs(w.diff))>(p*.001)){
	w.temp <- NR(w,x,z)$w.new
	w.diff <- w-w.temp
	w <- w.temp
	}
        w.new <- w
        lambda.new <- inv.logit(x %*% w.new)
	xbeta.new <- x %*% beta.new

        res <- (y - xbeta.new)^2

        s.new <- sqrt(sapply(1:k, function(i) sum(z[,i] * (res[, i]))/sum(z[, i])))

        sing <- sum(s.new < 1e-08)

        lambda <- lambda.new
        beta <- beta.new
        s <- s.new
        w <- w.new
        xbeta <- xbeta.new
        tmp <- t(-1/2/s[2]^2 * t(y - xbeta[, -1])^2) + 1/2/s[1]^2 *
            (y - xbeta[, 1])^2
        tmp2 <- (1-lambda)/s[2] * s[1]/lambda
        tmp3 <- log(1 + (tmp2 * exp(tmp)))
        newobsloglik <- sum(log(lambda)) - nlogroot2pi - n * log(s[1]) -
            1/2/s[1]^2 * sum((y - xbeta[, 1])^2) + sum(tmp3)


        if (sing > 0 || is.na(newobsloglik) || newobsloglik < 
            obsloglik || abs(newobsloglik) == Inf){# || sum(z) != n) {
            cat("Need new starting values due to singularity...", 
                "\n")
            restarts <- restarts + 1
            if (restarts > 15) 
                stop("Too many tries!")
	    tmp <- regmix.init(y = y, x = x, k = k, addintercept = addintercept)
	    beta = tmp$beta
	    s = tmp$s
	    k = tmp$k
	    w = rep(0,p*(k-1))
	    lambda = inv.logit(x %*% w)
	    diff <- 1
	    iter <- 0
	    xbeta <- x %*% beta
	    tmp <- t(-1/2/s[2]^2 * t(y - xbeta[, -1])^2) + 1/2/s[1]^2 *
	        (y - xbeta[, 1])^2
	    tmp2 <- (1-lambda)/s[2] * s[1]/lambda
	    tmp3 <- log(1 + (tmp2 * exp(tmp)))
	    obsloglik <- sum(log(lambda)) - nlogroot2pi - n * log(s[1]) -
	        1/2/s[1]^2 * sum((y - xbeta[, 1])^2) + sum(tmp3)
	    ll <- obsloglik
        }
	

	else{
        diff <- newobsloglik - obsloglik
        obsloglik <- newobsloglik
	ll <- c(ll, obsloglik)
        iter <- iter + 1
           if (verb) {
                cat("iteration=", iter, "diff=", diff, "log-likelihood", 
                  obsloglik, "\n")}
	}

    }
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
        rownames(beta) <- c(paste("beta", ".", 0:(p - 1), sep = ""))
        colnames(beta) <- c(paste("comp", ".", 1:k, sep = ""))
        colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
        a = list(x = x, y = y, w=w, lambda = cbind(lambda,1-lambda), beta = beta, 
            sigma = s, loglik = obsloglik, posterior = z, all.loglik = ll, 
            restarts = restarts, ft = "regmixEM")
        class(a) = "mixEM"
        a

}
