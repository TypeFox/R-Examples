# computing a likelihood value for the (E)CCC-GARCH(1,1) mdoel
   loglik.dcc1 <- function(param, dvar, model){
      nobs <- dim(dvar)[1]
      ndim <- dim(dvar)[2]
      In <- diag(ndim)
      param <- c(param, In[lower.tri(In)])
      para.mat <- p.mat(param, model, ndim)
      a <- para.mat$a
      A <- para.mat$A
      B <- para.mat$B
      h <- vector.garch(dvar, a, A, B)
      z <- dvar/sqrt(h)
      lf <- -0.5*nobs*ndim*log(2*pi)-0.5*sum(log(h))-0.5*sum(z^2)
      
      -lf
   }
