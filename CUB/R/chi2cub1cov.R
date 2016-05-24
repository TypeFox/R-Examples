#' @keywords internal


chi2cub1cov <-function(m,ordinal,covar,pai,gama){
  n<-length(ordinal)
  elle<-sort(unique(covar))
  
  kappa<-length(elle)
  
  matfrel<-matrix(NA,nrow=kappa,ncol=m)
  matprob<-matrix(NA,nrow=kappa,ncol=m)
  
  chi2<-0
  dev<-0
  
  j<-1
  while(j<=kappa){
    quali<-which(covar==elle[j])
    Wquali<-covar[quali]
    qualiord<-ordinal[quali]
    nk<- length(qualiord)
    matfrel[j,]=tabulate(qualiord,nbins=m)/nk
    nonzero<-which(matfrel[j,]!=0)
    paij<-pai
    csij<-1/(1+ exp(-gama[1]-gama[2]*elle[j]))
    matprob[j,]<-t(probcub00(m,paij,csij))
    chi2<-chi2+nk*sum(((matfrel[j,]-matprob[j,])^2)/matprob[j,])
    dev<- dev + 2*nk*sum(matfrel[j,nonzero]*log(matfrel[j,nonzero]/matprob[j,nonzero]))
    j<-j+1
  }
    
  df<-kappa*(m-1)-(length(gama)+1)
  cat("Degrees of freedom         ==>  df  =",df, "\n")
  cat("Pearson Fitting measure    ==>  X^2 =",chi2,"(p-val.=",1-pchisq(chi2,df),")","\n")
  cat("Deviance                   ==>  Dev =",dev,"(p-val.=",1-pchisq(dev,df),")","\n")
  results<-list('chi2'=chi2,'df'=df,'dev'=dev)
  
}
