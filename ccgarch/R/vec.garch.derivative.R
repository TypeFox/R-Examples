#************************************************************************************************************
   vec.garch.derivative <- function(dvar, B, h){     # the recursion for the derivative of h w.r.t. parameters
         nobs <- dim(h)[1]
         ndim <- dim(B)[1]
         npar.h <- ndim*(2*ndim+1)
         In <- diag(ndim)
         Iw <- diag(npar.h)                      # matrix of ones with order equal to the number of param in cond. var.
         s2 <- colMeans(dvar^2)
         v <- cbind(rep(1,nobs), dvar^2, h)

         vecIn <- as.vector(In)
         BIw = B%x%Iw                           # multiplier for the second term in the loop
         dhw <- matrix(0, ndim*npar.h, nobs)            # use columnwise
         dhw[,1] <- vecIn%x%c(1, rep(s2, 2))
         for (i in 2:nobs){
             dhw[,i] <- vecIn%x%v[(i-1),] + as.vector(BIw%*%dhw[,(i-1)])
         }
         t(dhw)                                    # (nobs x ndim*npar.h) matrix
   }
