bayesianLatentFit <- function(x,y,prob=0.01,c=100,iter=1000,burnin=NULL, siglevel = 0.05) {
    if(is.null(burnin)) burnin <- iter * 0.1

    p <- ncol(x)
    n <- nrow(x)

    GAMMA <- matrix(NA, ncol=p, nrow=burnin + iter)
    Z <- matrix(NA, ncol=n, nrow=burnin + iter)
    BETA <- matrix(NA, ncol=p, nrow=burnin + iter) 

    GAMMA[1,] <- rep(0,p)
    Z[1,] <- rep(0,n)
    BETA[1,] <- rep(0,p)

    for(i in 2:(burnin+iter)) {
        # Sample gamma_i values
        for(j in 1:p) {
            index <- which(GAMMA[i-1,] != 0 | (1:p) == j)
            xgamma <- x[,index]
            v1 <- exp(-0.5 * ( t(Z[i-1,]) %*% Z[i-1,] - (c)/(1+c) *
                   t(Z[i-1,]) %*% xgamma %*% solve(t(xgamma) %*% xgamma)  %*% t(xgamma) %*% Z[i-1,])) * prob
            #imat <- t(xgamma) %*% xgamma
            #diag(imat) <- diag(imat) + 1/c
            #v1 <- exp(-0.5 * ( t(Z[i-1,]) %*% Z[i-1,] - 
            #            t(Z[i-1,]) %*% xgamma %*% solve(imat)  %*% t(xgamma) %*% Z[i-1,])) * prob
            xgamma <- x[,setdiff(index,j)]
            if(length(xgamma) != 0) {
                v0 <- exp(-0.5 * ( t(Z[i-1,]) %*% Z[i-1,] - (c)/(1+c) *
                        t(Z[i-1,]) %*% xgamma %*% solve(t(xgamma) %*% xgamma)  %*% t(xgamma) %*% Z[i-1,])) * (1-prob)
                #imat <- t(xgamma) %*% xgamma
                #diag(imat) <- diag(imat) + 1/c
                #v0 <- exp(-0.5 * ( t(Z[i-1,]) %*% Z[i-1,] - 
                #        t(Z[i-1,]) %*% xgamma %*% solve(imat)  %*% t(xgamma) %*% Z[i-1,])) * (1-prob)
            } else {
                v0 <- exp(-0.5 * (t(Z[i-1,]) %*% Z[i-1,]) ) * (1-prob)
            }
            GAMMA[i,j] <- as.numeric( runif(1) <= v1 / (v1 + v0))
            if(is.na(GAMMA[i,j])) GAMMA[i,j] <- 0
        }

        # Sample z_i values
        mu <- x %*% BETA[i-1,]
        for(j in 1:n) {
            znew <- 0
            while(znew * sign((y[j] == 1) - 0.5) <= 0) {
                znew <- rnorm(1, mean = mu[j], sd=1)
            }
            Z[i,j] <- znew
        }

        # Sample beta_i values
        index <- which(GAMMA[i,] != 0 | (1:p) == j)
        if(length(index) != 0) {
            xgamma <- x[,index]
            V <- c / (1 + c) * solve(t(xgamma) %*% xgamma)
            BETA[i,-index] <- 0
            BETA[i,index] <- mvrnorm(n=1, mu=V %*% t(xgamma) %*% Z[i,], Sigma = V)
        } else {
            BETA[i,] <- 0
        }
    }

    return(list(beta = BETA, fitted = Z, model = GAMMA[-(1:burnin),]))

}
