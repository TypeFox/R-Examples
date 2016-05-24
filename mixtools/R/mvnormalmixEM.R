mvnormalmixEM = function (x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, arbmean = TRUE, arbvar = TRUE, 
    epsilon = 1e-08, maxit = 10000, verb = FALSE) 
{
    if(arbmean == FALSE && arbvar == FALSE){
    stop(paste("Must change constraints on mu and/or sigma!","\n"))
    }
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    tmp <- mvnormalmix.init(x = x, lambda = lambda, mu = mu, 
        sigma = sigma, k = k, arbmean=arbmean, arbvar = arbvar)
    lambda <- tmp$lambda
    mu<-tmp$mu
    sigma <- tmp$sigma
    k = tmp$k
    diff <- 1
    iter <- 0
    if (arbmean==FALSE){
        comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x, 
            mu, sigma[[i]]))
	} else{
    if (arbvar==FALSE) {
        comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x, 
            mu[[i]], sigma))
    }
    else comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x, 
        mu[[i]], sigma[[i]]))
    }
    comp <- sapply(comp, cbind)
    compsum <- apply(comp, 1, sum)
    obsloglik <- sum(log(compsum))
    ll <- obsloglik
	restarts <- 0
    while (diff > epsilon & iter < maxit) {
        if (arbvar) {
            z = matrix(nrow = n, ncol = k)
            for (i in 1:n) {
                for (j in 1:k) {
                  z.denom = c()
                  for (m in 1:k) {
                    z.denom = c(z.denom, lambda[m]/lambda[j] * 
                      (det(sigma[[j]])/det(sigma[[m]]))^(0.5) * 
                      exp(-0.5 * ((x[i, ] - mu[[m]]) %*% solve(sigma[[m]]) %*% 
                        t(t(x[i, ] - mu[[m]])) - (x[i, ] - mu[[j]]) %*% 
                        solve(sigma[[j]]) %*% t(t(x[i, ] - mu[[j]])))))
                  }
                  z[i, j] = 1/sum(z.denom)
                }
            }
	z = z/apply(z,1,sum)
#	  z[,k]=1-apply(as.matrix(z[,(1:(k-1))]),1,sum)
            sing <- sum(is.nan(z))
            lambda.new <- apply(z, 2, mean)
            if (sum(lambda.new < 1e-08)>0 || is.na(sum(lambda.new))) {
                sing <- 1
            }
            else {
            if(arbmean==FALSE) {
                  mu.new <- lapply(1:k, function(j) sapply(1:p, 
                  function(i) apply(z * x[, i], 2, sum))[j, ])
            mu.new <- apply(sapply(mu.new,as.vector),1,sum)/n
            mu.new <- lapply(1:k, function(j) mu.new)
            } else{
                mu.new <- lapply(1:k, function(j) sapply(1:p, 
                  function(i) apply(z * x[, i], 2, sum))[j, ]/sum(z[, 
                  j]))
            }
                sigma.new <- lapply(1:k, function(j) matrix(apply(sapply(1:n, 
                  function(i) z[i, j] * (x[i, ] - mu.new[[j]]) %*% 
                    t(x[i, ] - mu.new[[j]])), 1, sum), p, p)/sum(z[, 
                  j]))
                lambda <- lambda.new
                mu <- mu.new
                sigma <- sigma.new
                comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x, 
                  mu[[i]], sigma[[i]]))
                comp <- sapply(comp, cbind)
                compsum <- apply(comp, 1, sum)
                newobsloglik <- sum(log(compsum))
            }
        }
        else {
            z = matrix(nrow = n, ncol = k)
            sigma.inv = solve(sigma)
            for (i in 1:n) {
                for (j in 1:k) {
                  z.denom = c()
                  for (m in 1:k) {
                    z.denom = c(z.denom, lambda[m]/lambda[j] * 
                      (det(sigma.inv)/det(sigma.inv))^(0.5) * 
                      exp(-0.5 * ((x[i, ] - mu[[m]]) %*% sigma.inv %*% 
                        t(t(x[i, ] - mu[[m]])) - (x[i, ] - mu[[j]]) %*% 
                        sigma.inv %*% t(t(x[i, ] - mu[[j]])))))
                  }
                  z[i, j] = 1/sum(z.denom)
                }
            }
#	  z[,k]=1-apply(as.matrix(z[,(1:(k-1))]),1,sum)
	z = z/apply(z,1,sum)

            sing <- sum(is.nan(z))
            lambda.new <- apply(z, 2, mean)
            if (sum(lambda.new < 1e-08)>0 || is.na(sum(lambda.new))) {
                sing <- 1
            }
            else {
            if(arbmean==FALSE) {
                  mu.new <- lapply(1:k, function(j) sapply(1:p, 
                  function(i) apply(z * x[, i], 2, sum))[j, ])
            mu.new <- apply(sapply(mu.new,as.vector),1,sum)/n
            mu.new <- lapply(1:k, function(j) mu.new)
            } else{
                mu.new <- lapply(1:k, function(j) sapply(1:p, 
                  function(i) apply(z * x[, i], 2, sum))[j, ]/sum(z[, 
                  j]))
            }
                temp.sig <- lapply(1:k, function(j) matrix(apply(sapply(1:n, 
                  function(i) z[i, j] * (x[i, ] - mu.new[[j]]) %*% 
                    t(x[i, ] - mu.new[[j]])), 1, sum), p, p))
                sigma.new <- matrix(apply(sapply(temp.sig, as.vector), 
                  1, sum), p, p)/n
                lambda <- lambda.new
                mu <- mu.new
                sigma <- sigma.new
                comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x, 
                  mu[[i]], sigma))
                comp <- sapply(comp, cbind)
                compsum <- apply(comp, 1, sum)
                newobsloglik <- sum(log(compsum))
            }
        }
        if (sing > 0 || is.na(newobsloglik) || abs(newobsloglik) == Inf){# || sum(z) != n) {
            cat("Need new starting values due to singularity...", 
                "\n")
		restarts <- restarts + 1
		if(restarts>15) stop("Too many tries!")
            tmp <- mvnormalmix.init(x = x, k = k, arbmean=arbmean, arbvar = arbvar)
            lambda <- tmp$lambda
            mu <- tmp$mu
            sigma <- tmp$sigma
            k = tmp$k
            diff <- 1
            iter <- 0
            if (arbvar) {
                comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x, 
                  mu[[i]], sigma[[i]]))
            }
            else comp <- lapply(1:k, function(i) lambda[i] * 
                dmvnorm(x, mu[[i]], sigma))
            comp <- sapply(comp, cbind)
            compsum <- apply(comp, 1, sum)
            obsloglik <- sum(log(compsum))
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
    if(arbmean==FALSE) {
    mu = mu[[1]]
    }
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    cat("number of iterations=", iter, "\n")
    a=list(x=x, lambda = lambda, mu = mu, sigma = sigma, 
        loglik = obsloglik, posterior = z, all.loglik=ll, restarts=restarts, ft="mvnormalmixEM")
    class(a) = "mixEM"
    a
}
