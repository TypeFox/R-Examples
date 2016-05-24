# the objective function in the 2nd step DCC estimation. This is to be minimised!
loglik.dcc2 <- function(param, dvar){        # dvar is the standardised residuals
   nobs <- dim(dvar)[1]
   ndim <- dim(dvar)[2]
   DCC <- dcc.est(dvar, param)$DCC

#   lf <- numeric(ndim)  
   lf <- numeric(nobs)  # bug fixed on 2013.08.18
   for( i in 1:nobs){                        
      R <- matrix(DCC[i,], ndim, ndim)
      invR <- solve(R)
      lf[i] <- 0.5*(log(det(R)) +sum(dvar[i,]*crossprod(invR,dvar[i,])))  
   }
   sum(lf)
}






# the log-likelihood function for the 2nd step DCC estimation
#loglik.dcc2 <- function(param, dvar){        # dvar is the standardised residuals
#   nobs <- dim(dvar)[1]
#   ndim <- dim(dvar)[2]
#   
#   if(sum(param)>1|sum(param)<0){
#      param <- c(0, 0)
#   }
#   
#   DCC <- dcc.est(dvar, param)$DCC
#
#   lf <- numeric(ndim)
#   for( i in 1:nobs){                        
#      R <- matrix(DCC[i,], ndim, ndim)
#      invR <- solve(R)
#      lf[i] <- -0.5*(log(det(R)) +sum(dvar[i,]*crossprod(invR,dvar[i,])))
#   }
#   
#   sum(-lf)
#}





