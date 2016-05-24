summary.npmreg<-function(object,digits=max(3,getOption('digits')-3),...){
x<-object
k<-length(x$masses)
#Print formula
cat('\nCall: ',deparse(x$call),'\n\n')
#Print the estimated coefficients and SEs
cat('Coefficients:')
M1<-matrix(c(x$coefficients,x$SE.coefficients),ncol=2)
dimnames(M1)<-list(" "=c(names(x$coefficients))," "=c("Estimate ", "Std. Error"))
print.default(format(M1,digits=digits),print.gap=2,quote=FALSE,right=TRUE)
cat("\n")
#Print the estimated random effects distribution and SEs
if (k>1 & x$nrp==0){
  cat("Estimated NP Dist.:")
  M2<-matrix(c(x$mass.points,x$masses,x$SE.mass.points,x$SE.masses),ncol=2,nrow=(2*k))
  dimnames(M2)<-list(" "=c(paste('mass point',1:k,sep=' '),paste('mass',1:k,sep=' '))," "=c("Estimate", "Std. Error"))
  M3<-matrix(c(x$vcvremat,x$VRESE),nrow=1,ncol=2)
  dimnames(M3)<-list("Random effects variance: "=c(" ")," "=c("Estimate","Std. Error"))
  print.default(format(M2,digits=digits),print.gap=2,quote=FALSE,right=TRUE)
  cat("\n")
  cat("Random effects mean: constrained to 0","\n")
  print.default(format(M3,digits=digits),print.gap=2,quote=FALSE,right=TRUE)
}
if (k>1 & x$nrp>0){
  cat("Estimated Multivariate NP Dist.:","\n")
  cat("\n")
  cat("Mass Points:")
  M2<-matrix(c(x$mass.points,x$SE.mass.points),ncol=(1+x$nrp)*2,nrow=k)
  dimnames(M2)<-list(" "=c(paste('mass point',1:k,sep=' '))," "=c(dimnames(x$eBayes)[[2]],paste("Std. Error",dimnames(x$eBayes)[[2]])))
  print.default(format(M2,digits=digits),print.gap=2,quote=FALSE,right=TRUE)
  cat("\n")
  cat("Masses:")
  M3<-matrix(c(x$masses,x$SE.masses),ncol=2,nrow=k)
  dimnames(M3)<-list(" "=c(paste('mass',1:k,sep=' '))," "=c("Estimate", "Std. Error"))
  print.default(format(M3,digits=digits),print.gap=2,quote=FALSE,right=TRUE)
  cat("\n")
  cat("Random effects mean: all dimensions constrained to 0","\n","\n")
  cat("Random effects covariance (lower) and correlation (upper) matrices: ","\n")
  dimnames(x$var.cor.mat)<-list(" "=c(dimnames(x$eBayes)[[2]])," "=c(dimnames(x$eBayes)[[2]]))
  print.default(format(x$var.cor.mat,digits=digits),print.gap=2,quote=FALSE,right=TRUE)
  M4<-matrix(c(x$VRESE),nrow=1)
  dimnames(M4)<-list(paste("Std. Errors of random effects variances:")," "=dimnames(x$eBayes)[[2]])
  print.default(format(M4,digits=digits),print.gap=2,quote=FALSE,right=TRUE)
}
#Print -2*Log-likelihood
  cat("\n","-2(Log-Likelihood) :", format(signif(x$m2LogL,digits)),"\n")
  cat("Number of iterations :", x$iter,"\n")
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
