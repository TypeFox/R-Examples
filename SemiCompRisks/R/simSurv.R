
simSurv	<- function(cluster=NULL, x, beta.true, alpha.true, kappa.true, sigmaV.true=NULL, cens)
{
    if(!is.null(cluster) & is.null(sigmaV.true)){
        print("sigmaV.true must be given to simulate correated data")
    }
    else
    {
        n <- dim(x)[1]
        p <- dim(x)[2]
        
        if(is.null(cluster))
        {
            LP	<- as.vector(beta.true %*% t(x))
        }
        if(!is.null(cluster))
        {
            J <- length(unique(cluster))
            nj <- as.vector(table(cluster))
            
            V <- matrix(rnorm(J, 0, sqrt(sigmaV.true)), J, 1) # J X 1
            
            LP <- as.vector(beta.true %*% t(x) + rep(V[,1], nj))
        }
        
        T	<- rweibull(n, shape = alpha.true, scale = exp(-(log(kappa.true) + LP)/alpha.true))
        
        delta <- rep(NA, n)
        y		<- T
        Cen		<- runif(n, cens[1], cens[2])
        
        ind1	<- which(T < Cen)
        y[ind1] <- T[ind1]
        delta[ind1] <- 1
        
        ind0	<- which(T >= Cen)
        y[ind0] <- Cen[ind0]
        delta[ind0] <- 0
        
        ret <- data.frame(cbind(y, delta))
        return(ret)
    }
}
		

