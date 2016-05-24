spatialLmeUpdateMu <- function(I, G, V, cumV, mu0, betaBar, 
                     alpha,  sigma2inv, lambda2inv)
  .Call("spatialLmeUpdateMu", I, G, cumV, as.integer(sum(V)), mu0, betaBar, alpha,
        sigma2inv, lambda2inv)

spatialLmeUpdateSigma2invRate <- function(spatialMat, I, G, V, cumV, mu, alpha)
  .Call("spatialLmeUpdateSigma2invRate", spatialMat, I, G, cumV,
        as.integer(sum(V)), mu, alpha)

spatialLmeUpdateLambda2inv <- function(G, V, cumV, mu, mu0, c0, d0)
  .Call("spatialLmeUpdateLambda2inv", G, as.integer(V), cumV, mu, mu0, c0, d0)


spatial.lme.mcmc <- function(spatialMat, nlr, nsweep, verbose=TRUE){
    #-------------------------------------------------------------
    #Error checking
    if(nrow(spatialMat) != sum(nlr))
        stop("The number of locations within regions do not match the size of observations.")
    #-------------------------------------------------------------
    #Extract parameters from the input matrix

    #I: number of subjects within each condition
    I <- as.integer(ncol(spatialMat))
    G <- as.integer(length(nlr))
    cumVG <- as.integer(c(0, cumsum(nlr)))
     
    #Calculate values for some fixed value vectors
    #beta.bar.G: a vector of observed means of locations across subjects.
    #             The order is one region after another
    #
    beta.bar.G <- rowMeans(spatialMat)
    
    #beta.bar.I: a matrix of observed means of region within subjects.
    #            with one column per subject
    beta.bar.I <- apply(spatialMat, 2, function(x)
                        sapply(1:G, function(g) mean(x[(cumVG[g]+1):cumVG[g+1]])))
    #---------------------------------------------------------------
    #assign values to hyper-parameters
    #mu0.GJ: a prior of region means
    mu0.G <- sapply(1:G, function(g) mean(spatialMat[(cumVG[g]+1):cumVG[g+1],]))
    
    #set the hyper-parameters for wishart dist.
    h0 <- G
    
    #--------------------------------------------------------------
    #initialize mcmc results
    
    #lambda.inv.G: a matrix used to save precesions of region means
    lambda2.inv.G <- matrix(0, nrow=G, ncol=nsweep)
    lambda2.G  <- sapply(1:G, function(g) var(beta.bar.G[(cumVG[g]+1):cumVG[g+1]]))
    lambda2.inv.G[,1] <- 1 / lambda2.G
    
    #sigma2.inv.G: a matrix used to save precesions of residuals within each region
    sigma2.inv.G <- matrix(0, nrow=G, ncol=nsweep)
    sigma2.inv.G[,1] <- 1 / sapply(1:G, function(g)
                                   mean(sapply(1:I, function(i)
                                               var(spatialMat[(cumVG[g]+1):cumVG[g+1],i]))))
    
    #alpha.I: a matrix used to save regions effects per subject
    #within each subject, one region after another
    alpha.I <- matrix(0, nrow=I*G, ncol=nsweep)
    alpha.I[,1] <- as.vector(beta.bar.I) - rep(mu0.G, I)
   
    #mu.G: a matrix used to save location means 
    mu.G <- matrix(0, nrow=sum(nlr), ncol=nsweep)
    mu.G[,1] <-  beta.bar.G
    
    #Gamma.inv: an array used to save the inverse of variance covariance
    #matrix between regions
    #H0 is the hyper-parameter for wishart dist.
    Gamma.inv <- array(0, dim=c(G, G, nsweep))
    H0 <- matrix(0, G, G)
    H0 <- var(t(matrix( alpha.I[,1], nrow=G)))
    #Gamma.inv[,,1] <- chol2inv(chol(H0))
    Gamma.inv[,,1] <- solve(H0)
    #------------------------------------------------------------------------
    #assign values to the hyper-parameters a0, b0 for precisions of
    #of residuals within each region
    mu <- mean(sigma2.inv.G[,1])
    a0 <- 0.01
    b0 <- mu*100
    
    #assign values to the hyper-parameters c0, d0 of precesions of region means
    mu <- mean(lambda2.inv.G[,1])
    c0 <- 0.01
    d0 <- mu*100
 
    
    #--------------------------------------------------------------
    #run the simulation
    for(nsim in 2:nsweep){
        #update location means
        mu.G[,nsim] <- spatialLmeUpdateMu(I, G, nlr, cumVG, mu0.G, beta.bar.G,
                                alpha.I[,nsim-1],  
                                sigma2.inv.G[,nsim-1], lambda2.inv.G[,nsim-1])
          
        #update precisions for residuals within each region
        rate <- spatialLmeUpdateSigma2invRate(spatialMat, I, G, nlr, cumVG, 
                                    mu.G[,nsim], alpha.I[,nsim-1])
        rate <- b0 + 1/2*rate
        shape <- a0 + I*nlr/2
        sigma2.inv.g <- rgamma(G, shape=shape, rate=rate)
        sigma2.inv.G[,nsim] <- sigma2.inv.g
            
        #update regions effects per subject
        Gamma.inv.nsim <- Gamma.inv[,,nsim-1]
        D.inv <- diag((sigma2.inv.g * nlr))
        #Var <- chol2inv(chol(Gamma.inv.nsim + D.inv))
        Var <- solve(Gamma.inv.nsim + D.inv)
        mu.g <- mu.G[,nsim]
        mu.bar <- sapply(1:G, function(x) mean(mu.g[(cumVG[x]+1):cumVG[x+1]]))
        alpha.i <- NULL
        alpha.i.prod <- matrix(0, G, G)
        for(i in 1:I){
            Mu <- Var %*% (sigma2.inv.g * nlr * (beta.bar.I[,i] - mu.bar))
            alpha.i.i <- mvrnorm(1, Mu, Var)
            alpha.i <- c(alpha.i, alpha.i.i)
            alpha.i.prod <- alpha.i.prod +  alpha.i.i %*% t(alpha.i.i)
        }
        alpha.I[,nsim] <- alpha.i
                
        #update precesions of region means
        lambda2.inv.G[,nsim] <-  spatialLmeUpdateLambda2inv(G, nlr, cumVG,
                                                  mu.g, mu0.G, c0, d0)
            
        #update inverse of variance covariance matrix between regions
        #Gamma.inv[,,nsim] <- rwish(h0+I,  chol2inv(chol(h0*H0 + alpha.i.prod)))
        Gamma.inv[,,nsim] <- rwish(h0+I,  solve(h0*H0 + alpha.i.prod))
        
        if(verbose && nsim %% 1000 == 0){
            cat(paste(nsim, " iterations", " have finished.\n", sep=""))
        }
    }

    #Gamma <- sapply(1:nsim, function(i) chol2inv(chol(Gamma.inv[,,i])))
    Gamma <- sapply(1:nsim, function(i) solve(Gamma.inv[,,i]))

    list(mu.save=mu.G, sigma2.save=1/sigma2.inv.G, lambda2.save=1/lambda2.inv.G,
         Gamma.save=Gamma)
}
