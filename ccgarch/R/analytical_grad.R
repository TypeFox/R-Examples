# analytical score for ECCC
   analytical.grad <- function(a, A, B, R, u, model="diagonal"){
         invR <- solve(R)
      # defining number of length,dimensions, etc.      
         nobs <- dim(u)[1]            # number of observations
         ndim <- dim(u)[2]            # number of dimensions
         npar.h <- ndim*(2*ndim+1)      # number of param in GARCH eqns.
         npar.r <- ndim*(ndim-1)/2      # number of param in correlation matrix
         In <- diag(ndim)
         vecIn <- as.vector(In)
      # constructing volatilities
         h <- vector.garch(u, a, A, B)                # estimated volatility
            sq.invh <- 1/sqrt(h)
            norm.z <- u*sq.invh                         # standardized residuals
      # partial derivatives of D = diag(h^(1/2))
         dhw <- vec.garch.derivative(u,B,h)          # recursive equation for the partial derivatives of h w.r.t. parameters
         for (i in 1:ndim){
            dhw[,((i-1)*npar.h+1):(i*npar.h)] <- sq.invh[,i]*dhw[,((i-1)*npar.h+1):(i*npar.h)]
         }
         dhw <- t(0.5*dhw)                              # (ndim*npar.h x nobs) matrix

      if(model=="diagonal"){
         # ind <- c(rep(1,ndim), as.vector(In), as.vector(In))    This is wrong!!
         ind <- as.vector(rbind(1, In, In))
         dhw <- dhw[ind==1,]
         npar.h <- 3*ndim
      }
         vecdR <- vdR(ndim)                                        # partial derivatives of R w.r.t. correlation coefficients
         D <- matrix(0, npar.h, nobs)
         Q <- matrix(0, npar.r, nobs)
         for (i in 1: nobs){
            invD <- diag(sq.invh[i, ])
            Z <- norm.z[i, ]%o%norm.z[i, ]
            K <- Z%*%invR%*%invD
            d <- matrix(dhw[ ,i],npar.h,ndim)                       # npar.h x ndim
            W <- -0.5*as.vector(2*invD - (K + t(K)))                   # ndim^2 x 1
            W <- W[vecIn==1]                                       # removing rows multiplied by 0 in the next step. npar.h x ndim
            D[ ,i] <- d%*%W                                         # npar.h x nobs
            Q[ ,i] <- -0.5*vecdR%*%as.vector(invR - invR%*%Z%*%invR)  # npar.r x nobs
         }
      rbind(D,Q)
   }

