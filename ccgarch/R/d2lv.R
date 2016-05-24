d2lv <- function(u, B, h, model){ # the Hessian of the volatility part of the DCC log-likelihood function
      nobs <- dim(u)[1]
      ndim <- dim(u)[2]
      In <- diag(ndim)
      ind <- as.vector(rbind(1, In, In))
      dhw <- vec.garch.derivative(u, B, h)
      inv.h.sq <- 1/h^2
      if(model=="diagonal"){
         npar.h <- 3*ndim
      } else {
         npar.h <- ndim*(2*ndim+1)
      }
      dlv2 <- matrix(0, nobs, npar.h^2)
      for ( i in 1:nobs){
         dhwt <- matrix(dhw[i,],ncol=ndim)
         if(model=="diagonal"){
            dhwt <- dhwt[ind==1,]
         }
         tmp.dlv2 <- matrix(0, npar.h, npar.h)
         for ( j in 1:ndim){
           tmp.dlv2 <- tmp.dlv2 + inv.h.sq[i,j]*outer(dhwt[,j], dhwt[,j])
         }
         dlv2[i,] <- as.vector(tmp.dlv2)
      }
      0.5*dlv2
}
