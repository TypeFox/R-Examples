
simID <- function(cluster=NULL, x1, x2, x3, beta1.true, beta2.true, beta3.true,
                    alpha1.true, alpha2.true, alpha3.true,
                    kappa1.true, kappa2.true, kappa3.true, theta.true, SigmaV.true=NULL, cens)
{
    if(!is.null(cluster) & is.null(SigmaV.true)){
        print("SigmaV.true must be given to simulate correated data")
    }
    else
    {
        n <- dim(x1)[1]
        p1 <- dim(x1)[2]
        p2 <- dim(x2)[2]
        p3 <- dim(x3)[2]
        
        if(theta.true >0)
        {
            gamma.true <- rgamma(n, 1/theta.true, 1/theta.true)
        }
        if(theta.true == 0)
        {
            gamma.true <- rep(1, n)
        }
        
        
        if(is.null(cluster))
        {
            LP1	<- as.vector(beta1.true %*% t(x1))
            LP2	<- as.vector(beta2.true %*% t(x2))
            LP3	<- as.vector(beta3.true %*% t(x3))
        }
        if(!is.null(cluster))
        {
            J <- length(unique(cluster))
            nj <- as.vector(table(cluster))
            
            Vmat <- mvrnorm(J, rep(0, 3), SigmaV.true) # J X 3
            
            LP1 <- as.vector(beta1.true %*% t(x1) + rep(Vmat[,1], nj))
            LP2 <- as.vector(beta2.true %*% t(x2) + rep(Vmat[,2], nj))
            LP3 <- as.vector(beta3.true %*% t(x3) + rep(Vmat[,3], nj))
        }
        
        Rind <- NULL
        R <- rweibull(n, shape = alpha1.true, scale = exp(-(log(kappa1.true) +
        LP1 + log(gamma.true))/alpha1.true))
        D <- rweibull(n, shape = alpha2.true, scale = exp(-(log(kappa2.true) +
        LP2 + log(gamma.true))/alpha2.true))
        
        yesR <- R < D
        
        D[yesR] <- R[yesR] + rweibull(sum(yesR), shape = alpha3.true,
        scale = exp(-(log(kappa3.true) + LP3[yesR] + log(gamma.true[yesR]))/alpha3.true))
        delta1 <- rep(NA, n)
        delta2 <- rep(NA, n)
        y1 <- R
        y2 <- D
        Cen <- runif(n, cens[1], cens[2])
        ind01 <- which(D < R & D < Cen)
        y1[ind01] <- D[ind01]
        delta1[ind01] <- 0
        delta2[ind01] <- 1
        ind10 <- which(R < D & R < Cen & D >= Cen)
        y2[ind10] <- Cen[ind10]
        delta1[ind10] <- 1
        delta2[ind10] <- 0
        ind00 <- which(R >= Cen & D >= Cen)
        y1[ind00] <- Cen[ind00]
        y2[ind00] <- Cen[ind00]
        delta1[ind00] <- 0
        delta2[ind00] <- 0
        ind11 <- which(R < Cen & D < Cen & R < D)
        delta1[ind11] <- 1
        delta2[ind11] <- 1
        ret <- data.frame(cbind(y1, delta1, y2, delta2))
        return(ret)
    }
    
}


	