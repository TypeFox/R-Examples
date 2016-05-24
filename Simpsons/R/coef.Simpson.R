coef.Simpson <-
function(object, ...)
{
  Res <- data.frame(N = object$clustersize, Int = object$Allint,Beta = object$Allbeta, Pval=object$pvalues)
  group=cbind(object$totaln, object$groupint,object$groupbeta,object$pvaluesgr)
  rownames(group) <- ("Alldata")
  colnames(group)=colnames(Res)
  group=as.data.frame(group)
  rownames(Res) <- paste("Cluster",1:nrow(Res))
  Res=rbind(group,Res)
  return(Res)
}
