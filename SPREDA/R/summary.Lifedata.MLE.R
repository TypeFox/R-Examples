summary.Lifedata.MLE <-
function(object,...){
  m=length(object$coef)
  coef=c(object$coef[1:(m-1)], exp(object$coef[m]))
  names(coef)[m]="sigma"
  coef.dev=diag(m)
  coef.dev[m,m]=exp(object$coef[m])
  vcov=coef.dev%*%object$vcov%*%t(coef.dev)
  mat=matrix(0, nrow=m, ncol=4)
  mat[,1]=coef
  mat[,2]=sqrt(diag(vcov))
  mat[,3]=coef-1.96*mat[,2]
  mat[,4]=coef+1.96*mat[,2]
  mat[m,3]=exp(object$coef[m]-1.96*sqrt(diag(object$vcov)[m]))
  mat[m,4]=exp(object$coef[m]+1.96*sqrt(diag(object$vcov)[m]))
  colnames(mat)=c("mean", "std", "95% Lower", "95% Upper")
  rownames(mat)=names(coef)
  res=list(call=object$call, coef=coef, vcov=vcov, coefmat=mat, min=object$min, surv=object$surv, 
           dat=object$dat, ori.coef=object$coef, ori.vcov=object$vcov)  
  class(res)="summary.Lifedata.MLE"
  return(res)  
}
