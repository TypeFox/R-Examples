 sparsepredict=function(object,newx,s=NULL,type=c("response","coefficients","nonzero"),exact=FALSE,...){
 a0=t(as.matrix(object$a0))
 gamma=object$gamma
 gammarize=function(x,gamma){
   attr(x,"gamma")=gamma
   x
 }
 rownames(a0)="(Intercept)"
  nbeta=rbind2(a0,object$beta)
 
  if(!is.null(s)){
    vnames=dimnames(nbeta)[[1]]
    dimnames(nbeta)=list(NULL,NULL)
    lambda=object$lambda
    lamlist=lambda.interp(lambda,s)
    nbeta=nbeta[,lamlist$left,drop=FALSE]*lamlist$frac +nbeta[,lamlist$right,drop=FALSE]*(1-lamlist$frac)
    dimnames(nbeta)=list(vnames,paste(seq(along=s)))
  }
  if(type=="coefficients")return(gammarize(nbeta,gamma))
  if(type=="nonzero")return(gammarize(nonzeroCoef(nbeta[-1,,drop=FALSE],bystep=TRUE),gamma))
  gammarize(as.matrix(cbind2(1,newx)%*%nbeta),gamma)
  }
