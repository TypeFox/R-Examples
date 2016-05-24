###############################################################################
#######    r-values with nonparametric mixture dist.  #########################
###############################################################################


NPagrid <- function(estimate, nuisance, theta.alpha, theta.probs, 
                    alpha.grid, smooth, family, df = NULL) {
  
  ##########################################################################
  ##  Last edited: 1/5/14
  ##  Input
  ##    estimate - vector of local mles
  ##    nuisance - the associated vector of standard errors
  ##    support - the support of the estimated nonparametric mixture
  ##    mix_prop - the mixture proportions from the estimated nonparametric mixture 
  ##    smooth - the smoothing parameter for supsmu; used to smooth the
  ##                 empirical quantiles of Valpha
  ######
  ##   Output
  ##    a list with the following components:
  ##      rvalues - 
  ##      alpha - the alpha.grid vector (derived from mix_prop)
  ##      lam - The empirical quantiles of Valpha
  ##      lam_smooth - a smoothed version of lam
  ##      V - the "Valpha" matrix
  ###########################################################################
  
  nunits <- length(estimate)
  ngrid <- length(alpha.grid)
  
  rvalues <- mvec <- mmvec <- numeric(nunits)
  
  if(!is.null(df)) {
      if(length(df) == 1) {
          df <- rep(df, nunits)
      }
  }
  switch(family,
         gaussian={
             tmp <- PostProbNorm(estimate, nuisance, theta.alpha, theta.probs)
             V <- matrix(0, nrow = nunits, ncol = ngrid)
             for(i in 1:nunits) {
                  V[i,] <- cumsum(tmp$postprobs[i,])
             }
            # for(i in 1:nunits) {
            #      mvec[i] <- sum(dnorm(estimate[i],mean=theta.alpha,sd=nuisance[i])*theta.probs)
            # }
            # V <- matrix(0,nrow=nunits,ncol=ngrid)
            # V[,1] <- (theta.probs[1]/mvec)*(dnorm(estimate,mean=theta.alpha[1],sd=nuisance))
            # for(j in 2:ngrid)  {
            #     V[,j] <- V[,j-1] + (theta.probs[j]/mvec)*(dnorm(estimate,mean=theta.alpha[j],sd=nuisance))
            # }    
         },
         poisson={
             tmp <- PostProbPois(estimate, nuisance, theta.alpha, theta.probs)
             #PP <- tmp$postprobs
             V <- matrix(0,nrow=nunits,ncol=ngrid)
             ### theta.alpha goes from large to small
             for(i in 1:nunits) {
                 #mvec[i] <- sum(dpois(estimate[i],lambda=nuisance[i]*theta.alpha)*theta.probs)
                 V[i,] <- cumsum(tmp$postprobs[i,])
             }

         },
         binomial={
             tmp <- PostProbBinomial(estimate, nuisance, theta.alpha, theta.probs)
             V <- matrix(0, nrow = nunits, ncol = ngrid)
             for(i in 1:nunits) {
                  V[i,] <- cumsum(tmp$postprobs[i,])
             }
         },
         tdist={
              tmp <- PostProbT(estimate, nuisance, df, theta.alpha, theta.probs)
              V <- matrix(0, nrow = nunits, ncol = ngrid)
              for(i in 1:nunits) {
                   V[i,] <- cumsum(tmp$postprobs[i,])
              }      
              #for(i in 1:nunits) {
              #    mvec[i] <- sum((dt((estimate[i] - theta.alpha)/nuisance[i],df=df[i])/nuisance[i])*theta.probs)
              #}
              #V <- matrix(0,nrow=nunits,ncol=ngrid)
              #V[,1] <- (theta.probs[1]/mvec)*(dt((estimate - theta.alpha[1])/nuisance,df=df)/nuisance)
              #for(j in 2:ngrid)  {
              #    V[,j] <- V[,j-1] + (theta.probs[j]/mvec)*(dt((estimate - theta.alpha[j])/nuisance,df=df)/nuisance)
              #}
         },
   )
  
  lam <- numeric(ngrid)
  for( j in 1:ngrid) {
     lam[j] <- quantile(V[,j], prob= 1 - alpha.grid[j], names = FALSE, type = 1)
  }

  if(smooth=="none")  {
     #cc2 <- approxfun(alpha.grid,lam,yleft=1,yright=0)
     lam.smooth.eval <- lam
     lam.smooth <- approxfun( c(0,alpha.grid,1), c(1,lam,0))
  }
  else {
     cc2 <- supsmu(alpha.grid, lam, bass = smooth)
     lam.smooth.eval <- cc2$y
     lam.smooth <- approxfun( c(0,cc2$x,1), c(1,cc2$y,0))
  }
      
  ### For each row of Valpha, determine the index at which
  ### Valpha[i,] intersects lambda_{\alpha}
  #cut.ind <- VVcut(V, lam.smooth.eval, nrow(V), length(lam), alpha.grid) 
  rvalues <- VVcut(V, lam.smooth.eval, nrow(V), length(lam), alpha.grid) 
  ## Interpolation:
  #rvalues <- rep(0,nunits)
  #for(j in 1:nunits) {
     # print(cut.ind[j])
  #    if(cut.ind[j] !=1 ) {
  #         g0 = V[j,cut.ind[j] - 1] - lam.smooth.eval[cut.ind[j]-1] 
  #         gdelt = V[j,cut.ind[j]] - lam.smooth.eval[cut.ind[j]] - g0
  #         slop = (alpha.grid[cut.ind[j]] - alpha.grid[cut.ind[j] - 1])/gdelt
  #         rvalues[j] = alpha.grid[cut.ind[j] - 1] - g0*slop
  #    }
  #    else {
  #        rvalues[j] <- alpha.grid[1]
  #     }
  #} 
  #rvalues <- alpha.grid[cut.ind]
  #rvalues <- cut.ind

  ans <- list()
  ans$rvalues <- rvalues
  ans$alpha <- alpha.grid
  ans$lam <- lam
  ans$lamsmooth <- lam.smooth
  ans$V <- V
  return(ans)
}


