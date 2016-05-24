poisgam <- function(X, y, start.coef,
                    start.sigma, fix.sigma,
                    cluster, offset){

    Pfun <- function(par, X, y, cluster, offset = NULL){
        ## Poisson-gamma; likelihood function
        ## Last component in 'par' is the gamma -log(variance)
        ## First component in par and beta is constant == 1 (Intercept).
        
        if (is.null(offset)) offset <- rep(0, length(y))
        
        p <- length(par) - 1
        alpha <- par[p + 1]
        beta <- par[-(p + 1)]
        ealpha <- exp(alpha)
        cls <- unique(cluster)
        n <- length(cls)
        score <- offset + X %*% beta
        escore <- exp(score)
        fun <- sum(score * y) ## Can be made more efficient (later)
        if (alpha > 600){ # "ZERO" variance
            cat("ZERO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
            fun <- fun - sum(escore) - sum(lgamma(y + 1))
        }else{
            for (i in cls){
                who <- (cluster == i)
                sumy <- sum(y[who])
                sums <- sum(escore[who])
                ##cat("alpha = ", alpha, "\n")
                ##cat("ealpha = ", ealpha, "\n")
                ##cat("sumy = ", sumy, "\n")
                ##cat("i = ", i, "\n\n")
                fun <- fun + lgamma(ealpha + sumy) + alpha * ealpha -
                    lgamma(ealpha) - (ealpha + sumy) * log(ealpha + sums) -
                        sum(lgamma(y[who] + 1))
            }
        }
        fun
    }

    Pgrad <- function(par, X, y, cluster, offset = NULL){
        ## Poisson-gamma; gradient function
        ## Last component in 'par' is the gamma -log(variance)
        ## First component in par and beta is constant == 1 (Intercept).
        ## Returns gradient vector (p + 1)
        
        if (is.null(offset)) offset <- rep(0, length(y))
        p <- length(par) - 1
        alpha <- par[p + 1]
        beta <- par[-(p + 1)]
        ealpha <- exp(alpha)
        cls <- unique(cluster)
        n <- length(cls)
        score <- offset + X %*% beta
        grad <- numeric(p + 1)
        grad[p + 1] <- n * (alpha + 1 - digamma(ealpha))
        for (m in 1:p){
            grad[m] <- sum(y * X[, m])
        }
        escore <- exp(score)
        for (i in cls){
            who <- (cluster == i)
            sumy <- sum(y[who])
            sums <- sum(escore[who])
            for (m in 1:p){
                grad[m] <- grad[m] -
                    (ealpha + sumy) * sum(X[who, m] * escore[who]) /
                        (ealpha + sums)
            }
            grad[p + 1] <- grad[p + 1] + digamma(ealpha + sumy) -
                log(ealpha + sums) - (ealpha + sumy) / (ealpha + sums)
        }
        grad[p + 1] <- ealpha * grad[p + 1] # Note!
        grad
    }

    p <- NCOL(X)
    
    par <- rep(0, p + 1)
    cat("res0: \n")
    res0 <- optim(par, Pfun, gr = NULL,
                  X = X, y = y, cluster = cluster, offset = offset,
                  ##method = "L-BFGS-B",
                  ##method = "BFGS",
                  ##lower = rep(-50, p + 1),
                  ##upper = rep(50, p + 1),
                  control = list(fnscale = -1, trace = 5), hessian = TRUE)
    cat("res1: \n")
    res1 <- optim(par, Pfun, gr = Pgrad,
                  X = X, y = y, cluster = cluster, offset = offset,
                  ##method = "L-BFGS-B",
                  method = "BFGS",
                  ##lower = rep(-50, p + 1),
                  ##upper = rep(50, p + 1),
                  control = list(fnscale = -1, trace = 5), hessian = TRUE)
    
    vari <- solve(-res1$hessian)
    alpha <- res1$par[p + 1]
    sigma <- exp(-alpha)
    sds <- sqrt(diag(vari))
    sigma.sd <- sqrt(sds[p + 1]) / alpha
    beta <- res1$par[-(p + 1)]
    list(beta = beta,
         beta.sd = sds[-(p + 1)],
         sigma = sigma,
         sigma.sd = sigma.sd,
         loglik = res1$value,
         variance = vari,
         bootP = NULL,
         post.mode = NULL,
         post.mean = NULL
         ) 
}


    
        
