# Analytical Hessian in ECCC-GARCH
   analytical.Hessian <- function(a, A, B, R, u, model){
      nobs <- dim(u)[1]              # number of observations
      ndim <- dim(u)[2]              # number of dimensions
      npar.h <- ndim*(2*ndim+1)      # number of param in GARCH eqns.
      npar.r <- ndim*(ndim-1)/2      # number of param in correlation matrix
      In <- diag(ndim)
      vecIn <- as.vector(In)
      invR <- solve(R)
   # constructing volatilities
      h <- vector.garch(u, a, A, B)
         sq.h <- sqrt(h)
         sq.invh <- 1/sq.h
   # partial derivatives of D = diag(h^(1/2))
      dhw <- vec.garch.derivative(u,B,h)                     # recursive equation for the partial derivatives of h w.r.t. parameters
                                                   # dim(dhw) = (nobs, ndim*npar.h)
      for (i in 1:ndim){
         dhw[,((i-1)*npar.h+1):(i*npar.h)] <- sq.invh[,i]*dhw[,((i-1)*npar.h+1):(i*npar.h)]
      }
      dhw <- t(0.5*dhw)                              # (ndim*npar.h x nobs) matrix 

      if(model=="diagonal"){
         ind <- c(rep(1,ndim), as.vector(In), as.vector(In))
         dhw <- dhw[ind==1,]
         npar.h <- 3*ndim
      }
         vecdR <- vdR(ndim)                                        # partial derivatives of R w.r.t. correlation coefficients
         H11 <- matrix(0,npar.h,npar.h)
         H21 <- matrix(0,npar.r,npar.h)
         for (i in 1:nobs){
            D <- diag(sq.h[i,])
            invD <- diag(sq.invh[i,])
            invH <- invD%*%invR%*%invD
            dh <- matrix(dhw[,i],npar.h,ndim)                      # dvec(D_t)'/dw
            H11.temp <- 2*(invD%x%invD) + invH%x%R + R%x%invH
            H11.temp <- H11.temp[vecIn==1,vecIn==1]                # removing rows and columns multiplied by 0 in the next step. npar.h x ndim
            H11 <- H11 + dh%*%H11.temp%*%t(dh)
            RD <- invR%*%invD
            H21.temp <- vecdR%*%(RD%x%In + In%x%RD)
            H21.temp <- H21.temp[,vecIn==1]
            H21 <- H21 + H21.temp%*%t(dh)
         }

         H22 <- nobs*vecdR%*%(invR%x%invR)%*%t(vecdR)
         0.5*cbind(rbind(H11,H21),rbind(t(H21),H22))         # analytical Hessian
   }
