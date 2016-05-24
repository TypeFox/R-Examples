# the full log-likelihood function of the DCC-GARCH(1,1) model
loglik.dcc <- function (param, dvar, model) {
    nobs <- dim(dvar)[1]
    ndim <- dim(dvar)[2]
    Id <- diag(ndim)
    npar <- length(param)
    dcc.param <- param[(npar-1):npar]

      if(model=="diagonal"){                       # for the diagonal vector GARCH equation
         a <- param[1:ndim]                         # constant in variance
         A <- diag(param[(ndim+1):(2*ndim)])        # ARCH parameter
         B <- diag(param[(2*ndim+1):(3*ndim)])      # GARCH parameter
      } else if(model=="extended"){                # for the extended vector GARCH equation
         a <- param[1:ndim]
         A <- matrix(param[(ndim+1):(ndim^2+ndim)], ndim, ndim)
         B <- matrix(param[(ndim^2+ndim+1):(2*ndim^2+ndim)], ndim, ndim)
      } 

    h <- vector.garch(dvar, a, A, B)
    z <- dvar/sqrt(h)
    DCC <- dcc.est(z, dcc.param)$DCC
    
#    lf1 <- -0.5*ndim*log(2*pi) - rowSums(log(h)) -0.5*rowSums(z^2)
    lf1 <- -0.5*ndim*log(2*pi) - rowSums(log(h))         # bug fixed on 23072009
    lf2 <- numeric(nobs)
   for( i in 1:nobs){                        
      R <- matrix(DCC[i,], ndim, ndim)
      invR <- solve(R)
      tmpz <- as.vector(z[i,])
      lf2[i] <- -0.5*(log(det(R)) +sum(tmpz*crossprod(invR, tmpz)))
   }
    sum(lf1 + lf2)
}
