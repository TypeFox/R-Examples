predict.sparsenet=function(object,newx,s=NULL,which.gamma=NULL,type=c("response","coefficients","nonzero"),exact=FALSE,...){
 type=match.arg(type)
  if(missing(newx)){
    if(!match(type,c("coefficients","nonzero"),FALSE))stop("You need to supply a value for 'newx'")
     }
  coeflist=object$coefficients
  ngamma=length(coeflist)
   coeflistseq=seq(along=coeflist)
 if(is.null(which.gamma))which.gamma=coeflistseq
 else   which.gamma=coeflistseq[match(which.gamma,coeflistseq,0)]
if(length(which.gamma)>1){
  predlist=as.list(which.gamma)
  names(predlist)=names(coeflist)[which.gamma]
  for(j in seq(along=which.gamma))predlist[[j]]=sparsepredict(coeflist[[which.gamma[j]]],newx,s,type,...)
  predlist
}
 else sparsepredict(coeflist[[which.gamma]],newx,s,type,...)
}
