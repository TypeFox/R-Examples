# numerical gradient of the full log-likelihood of the (E)DCC-GARCH model.
grad.dcc.full <- function(a, A, B, dcc.para, dvar, d=1e-5, model){
   if(model=="diagonal"){
      param <- c(a, diag(A), diag(B), dcc.para)
   }  else  {
      param <- c(a, as.vector(A), as.vector(B), dcc.para)
   }
   nobs <- dim(dvar)[1]
   ndim <- dim(dvar)[2]
   npara <- length(param)
   Id <- d*diag(npara)
   
   param.d <- matrix(param, npara, npara)+Id
   
   lf <- loglik.dcc(param, dvar, model)
   lf.d <- matrix(0, nobs, npara)
   for(i in 1:npara){
      lf.d[,i] <- loglik.dcc(param.d[,i], dvar, model)
   }
   (lf.d - lf)/d
}
