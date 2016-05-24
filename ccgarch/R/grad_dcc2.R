# numerical gradient of the objective function in the 2nd step DCC estimation.
grad.dcc2 <- function(param, dvar, d=1e-5){
   nobs <- dim(dvar)[1]
   ndim <- dim(dvar)[2]
   npara <- length(param)
   Id <- d*diag(npara)
   param1 <- param + Id[,1]
   param2 <- param + Id[,2]

#   for(i in 1:npara){
#      assign(paste("param", i, sep=""), param+Id[,i])
#   }
   DCC <- dcc.est(dvar, param)$DCC
   DCC1 <- dcc.est(dvar, param1)$DCC
   DCC2 <- dcc.est(dvar, param2)$DCC

   lf <- numeric(ndim)
   lf1 <- numeric(ndim)
   lf2 <- numeric(ndim)
   for( i in 1:nobs){                        
      R <- matrix(DCC[i,], ndim, ndim)
      R1 <- matrix(DCC1[i,], ndim, ndim)
      R2 <- matrix(DCC2[i,], ndim, ndim)
      invR <- solve(R)
      invR1 <- solve(R1)
      invR2 <- solve(R2)
      lf[i] <- 0.5*(log(det(R)) +sum(dvar[i,]*crossprod(invR, dvar[i,])))
      lf1[i] <- 0.5*(log(det(R1)) +sum(dvar[i,]*crossprod(invR1, dvar[i,])))
      lf2[i] <- 0.5*(log(det(R2)) +sum(dvar[i,]*crossprod(invR2, dvar[i,])))
   }
   c(sum((lf1 - lf)/d), sum((lf2 - lf)/d) )
}
