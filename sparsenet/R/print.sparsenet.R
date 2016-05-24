print.sparsenet=function(x,digits = max(3, getOption("digits") - 3),...){
  cat("\nCall: ", deparse(x$call), "\n\n")
  coeflist=x$coefficients
  rsq=x$rsq
  gnames=names(coeflist)
  for(i in seq(along=coeflist)){
    x=coeflist[[i]]
    cat(gnames[i],": Gamma = ",signif(x$gamma,digits),"\n")
    print(cbind(Df=x$df,"Rsq"=signif(rsq[i,],digits),Lambda=signif(x$lambda,digits)))
  }
}
