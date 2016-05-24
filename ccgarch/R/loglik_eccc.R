# computing a likelihood value for the (E)CCC-GARCH(1,1) mdoel
   loglik.eccc <- function(param, dvar, model){
      nobs <- dim(dvar)[1]
      ndim <- dim(dvar)[2]
      para.mat <- p.mat(param, model, ndim)
      a <- para.mat$a
      A <- para.mat$A
      B <- para.mat$B
      R <- para.mat$R
      
      # check if R is positive definite
      eigenR <- eigen(R)$values
      if(max(abs(R[lower.tri(R)]))>1.0||min(eigenR)<0||!is.double(eigenR)){
         R <- diag(ndim)
      }
      h <- vector.garch(dvar, a, A, B)
      z <- dvar/sqrt(h)
      lndetR <- log(det(R))
      invR <- solve(R)
      lf <- -0.5*nobs*ndim*log(2*pi) - 0.5*sum(log(h)) - 0.5*nobs*lndetR - 0.5*sum((z%*%invR)*z)
      -lf
   }
