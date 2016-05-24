# constructing parameter vector and matrices. positivity constraints are imposed in "lleccc.c"
   p.mat <- function(para, model, ndim){
         npara <- length(para)
      if(model=="diagonal"){                         # for the diagonal vector GARCH equation
         a <- para[1:ndim]^2                          # constant in variance
         A <- diag(para[(ndim+1):(2*ndim)]^2)               # ARCH parameter
         B <- diag(para[(2*ndim+1):(3*ndim)]^2)                   # GARCH parameter
         R <- diag(ndim)                              # Constant Conditional Correlation Matrix
         R[lower.tri(R)] <- para[(3*ndim+1):npara]; R <- (R+t(R)); diag(R) <- 0.5*diag(R)
      } else if(model=="extended"){   # for the extended vector GARCH equation
         a <- para[1:ndim]^2
         A <- matrix(para[(ndim+1):(ndim^2+ndim)]^2, ndim, ndim)
         B <- matrix(para[(ndim^2+ndim+1):(2*ndim^2+ndim)]^2, ndim, ndim)
         R <- diag(ndim)
         R[lower.tri(R)] <- para[(2*ndim^2+ndim+1):npara]; R <- (R+t(R)); diag(R) <- 0.5*diag(R)
      } else if(model=="ECCC.neg"){   # the extended model with negative interaction
         a <- para[1:ndim]^2
         A <- matrix(para[(ndim+1):(ndim^2+ndim)]^2, ndim, ndim)
         B <- matrix(para[(ndim^2+ndim+1):(2*ndim^2+ndim)], ndim, ndim)
         R <- diag(ndim)
         R[lower.tri(R)] <- para[(2*ndim^2+ndim+1):npara]; R <- (R+t(R)); diag(R) <- 0.5*diag(R)
      }
      list(a=a, A=A, B=B, R=R)
   }
