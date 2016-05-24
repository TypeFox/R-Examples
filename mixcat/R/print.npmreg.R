print.npmreg<-function(x,digits=max(3,getOption('digits')-3), ...){
   k<-length(x$masses)
   #Print formula
   cat('\nCall: ',deparse(x$call),'\n\n')
   #Print the estimated coefficients
   cat('Coefficients:\n')
   print.default(format(x$coefficients,digits=digits),print.gap=2,quote=FALSE)
   cat("\n")
   #Print the estimated random effects distribution
   if (k>1 & x$nrp==0){
      cat("Estimated NP Dist.:\n")
      print.default(format(c(x$mass.points,x$masses),digits=digits),print.gap=2,quote=FALSE)
      cat("\n")
      cat("Random intercept mean: constrained to 0","\n")
      cat("Random intercept variance: ",format(signif(x$var.cor.mat,digits)),"\n")}
   if (k>1 & x$nrp>0){
      cat("Estimated Multivariate NP Dist.:")
      M2<-matrix(c(x$mass.points,x$masses),ncol=(x$nrp+2),nrow=k)
      dimnames(M2)<-list(" "=c(paste('mass point',1:k,sep=' '))," "=c(dimnames(x$eBayes)[[2]],"mass"))
      print.default(format(M2,digits=digits),print.gap=2,quote=FALSE,right=TRUE)
      cat("\n")
      cat("Random effects mean : all dimensions constrained to 0","\n","\n")
      dimnames(x$var.cor.mat)<-list(" "=c(dimnames(x$eBayes)[[2]])," "=c(dimnames(x$eBayes)[[2]]))
      cat("Random effects covariance (lower) and correlation (upper) matrices: ","\n")
      print.default(format(x$var.cor.mat,digits=digits),print.gap=2,quote=FALSE,right=TRUE)
      }
   #Print -2*Log-likelihood
   cat("\n","-2(Log-Likelihood) :", format(signif(x$m2LogL,digits)),"\n")
   if (x$iter > (x$maxit-1) | x$flagcvm == x$iter | x$flaginfo == x$iter) cat("\n")
   #Conditional print of warning about number of iterations
   if (x$iter > (x$maxit-1))
      cat("Warning: maximum mumber of iterations was reached","\n")
   #Conditional print of flag cvm
   if (x$flagcvm == x$iter)
      cat("Warning: eigenvalues of multinomial covariance matrix less than argument 'tol' appear at last  iteration","\n")
   #Conditional print of flag info
   if (x$flaginfo == x$iter)
      cat("Warning: eigenvalues of information matrix less than argument 'tol' appear at last iteration", "\n")
}
