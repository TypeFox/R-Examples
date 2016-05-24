repnormmixEM = function (x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, arbmean = TRUE, arbvar = TRUE, 
    epsilon = 1e-08, maxit = 10000, verb = FALSE) 
{
    if(arbmean == FALSE && arbvar == FALSE){
   stop(paste("Must change constraints on mu and/or sigma!","\n"))
    }
    s = sigma
    n <- ncol(x)
    m <- nrow(x)
    tmp <- repnormmix.init(x = x, lambda = lambda, mu = mu, s = s, 
        k = k, arbmean = arbmean, arbvar = arbvar)
    lambda <- tmp$lambda
    mu <- tmp$mu
    s <- tmp$s
    s.2 <- s^2
    k <- tmp$k
    diff <- 1
    iter <- 0
    res <- lapply(1:k, function(i) (x - mu[i * arbmean + (1 - arbmean)])^2)
    comp <- lapply(1:k, function(i) lambda[i] * (2 * pi * s.2[i * 
        arbvar + (1 - arbvar)])^(-m/2) * exp(-0.5 * apply(res[[i]], 
        2, sum)/s.2[i * arbvar + (1 - arbvar)]))
    obsloglik <- sum(log(apply(sapply(comp, rbind), 1, sum)))
    ll <- obsloglik
	restarts <- 0
    while (diff > epsilon & iter < maxit) {
        xmu <- lapply(1:k, function(i) apply(res[[i]], 2, sum))
        z = matrix(nrow = n, ncol = k)
        for (i in 1:n) {
            for (j in 1:k) {
                z.denom = c()
                for (h in 1:k) {
                  z.denom = c(z.denom, lambda[h]/lambda[j] * 
                    (s[j * arbvar + (1 - arbvar)]/s[h * arbvar + 
                      (1 - arbvar)])^m * exp(-0.5 * ((1/s[h * 
                    arbvar + (1 - arbvar)]^2) * xmu[[h]][i] - 
                    (1/s[j * arbvar + (1 - arbvar)]^2) * xmu[[j]][i])))
                }
                z[i, j] = 1/sum(z.denom)
            }
        }
	z = z/apply(z,1,sum)
#	  z[,k]=1-apply(as.matrix(z[,(1:(k-1))]),1,sum)

        lambda.new <- apply(z, 2, mean)
        if (sum(lambda.new < 1e-08)>0 || is.na(sum(lambda.new))) {
            sing <- 1
        }
        else {
        if(arbmean){
            mu.new <- sapply(1:k, function(i) sum(t(z[, i] * 
                t(x))))/(m * apply(z, 2, sum))
        } else {
            mu.new <- sum(sapply(1:k, function(i) sum(t(z[, i] * 
                t(x)))))/(m * n)
        }
            res <- lapply(1:k, function(i) (x - mu.new[i * arbmean + (1 - arbmean)])^2)
            if (arbvar) {
                s.new <- sqrt(sapply(1:k, function(i) sum(t(z[, 
                  i] * t(res[[i]]))))/(m * apply(z, 2, sum)))
            }
            else s.new <- sqrt(sum(sapply(1:k, function(i) sum(t(z[, 
                i] * t(res[[i]])))))/(m * n))
            lambda <- lambda.new
            mu <- mu.new
            s <- s.new
            s.2 <- s^2
            sing <- sum(is.nan(z))
            comp <- lapply(1:k, function(i) lambda[i] * (2 * 
                pi * s.2[i * arbvar + (1 - arbvar)])^(-m/2) * 
                exp(-0.5 * apply(res[[i]], 2, sum)/s.2[i * arbvar + 
                  (1 - arbvar)]))
            newobsloglik <- sum(log(apply(sapply(comp, rbind), 
                1, sum)))
        }
        if (sing > 0 || is.na(newobsloglik) || abs(newobsloglik) == Inf ){#|| sum(z) != n) {
            cat("Need new starting values due to singularity...", 
                "\n")
		restarts <- restarts + 1
		if(restarts>15) stop("Too many tries!")
            tmp <- repnormmix.init(x = x, k = k, arbmean = arbmean, arbvar = arbvar)
            lambda <- tmp$lambda
            mu <- tmp$mu
            s <- tmp$s
            s.2 <- s^2
            k <- tmp$k
            diff <- 1
            iter <- 0
            res <- lapply(1:k, function(i) (x - mu[i * arbmean + (1 - arbmean)])^2)
            comp <- lapply(1:k, function(i) lambda[i] * (2 * 
                pi * s.2[i * arbvar + (1 - arbvar)])^(-m/2) * 
                exp(-0.5 * apply(res[[i]], 2, sum)/s.2[i * arbvar + 
                  (1 - arbvar)]))
            obsloglik <- sum(log(apply(sapply(comp, rbind), 1, 
                sum)))
        ll <- obsloglik
        }
        else {
            diff <- newobsloglik - obsloglik
            obsloglik <- newobsloglik
        ll <- c(ll, obsloglik)
            iter <- iter + 1
            if (verb) {
                cat("iteration=", iter, "diff=", diff, "log-likelihood", 
                  obsloglik, "\n")
            }
        }
    }
    scale.order = order(s)
    sigma.min = min(s)
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    if(arbmean == FALSE){
    z = z[,scale.order]
colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=as.data.frame(x), lambda = lambda[scale.order], mu = mu, sigma = sigma.min, scale = s[scale.order]/sigma.min, loglik = obsloglik, 
        posterior = z[,scale.order], all.loglik=ll, restarts = restarts, ft="repnormmixEM")
    class(a) = "mixEM"
    a
    } else {
colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=as.data.frame(x), lambda = lambda, mu = mu, sigma = s, loglik = obsloglik, 
        posterior = z, all.loglik=ll, restarts = restarts, ft="repnormmixEM")
    class(a) = "mixEM"
    a
    }
}
