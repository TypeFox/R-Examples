"predict.spa" <-function(object,xnew,gnew,type=c("vector","probs","both"),...){
  if(class(object)!="spa"){
    stop("Error: object must be of type spa")
  }
  if(missing(type)){
    type="vector"
}
  if(missing(gnew))
    return(fitted(object))
  
  if(!missing(gnew)){
    if(is.data.frame(gnew))
      gnew=as.matrix(gnew)
    if(!is.matrix(gnew))
      stop(paste("Error in graph:  must be of type  matrix"))
  }
  p=0
  if(!missing(xnew)){
    xnew=as.matrix(xnew)
    if(dim(xnew)[1]!=n)
      stop("Error in dims x: dim(x)[1]!=dim(gnew)[1]")
    p=dim(xnew)[2]
    if(dim(object)[3]!=p){
      stop("Error:  The dims x do not match: Expect: (",dim(object)[3],") Inputted: (",p,")")
    }
  }
  n=dim(gnew)[1]
  np=dim(gnew)[2]
  if(np<1){
    stop("Error: graph must be a matrix\n")
  }
  if(np>1){
    if(dim(object)[1]!=np){
      stop("Error:  The dims graph do not match: Expect: (",dim(object)[1],") Inputted: (",np,")")
    }
  }
  xdat=!is.null(object$model$xstr)
  if(missing(xnew)& xdat){
    stop("Error: object x must be supplied whenever object was fit with x data")
  }
  ctype=object$type
  yvec<-fitted(object)
  yvec[object$model$L]=object$model$y
  W=gnew
  if(object$control$dissimilar)
    W=object$kern(gnew,object$model$parm.est$cvlam)+object$control$adjust
  if(!xdat){
    if(np==1){
      fit=sum(W*yvec)/sum(W)
    }else{
      S=W/apply(W,1,sum)
      fit=as.vector(S%*%yvec)
    }
  }else{
    if(object$control$warn)
      warning("Inductive prediction of semi-par graph-based models is still under development")
    bet=coef(object)
    S= S-t(matrix(apply(S,2,sum),np,n))/n
    f2=S%*%yvec-(S%*%xnew)%*%bet
    f1=xnew%*%bet
    fit=as.vector(f2+f1)
  }
  if(ctype=="soft"){
    return(fit)
  }
  ty<-as.factor(fit>0.5)
  if(type=="vector"){
    return(ty)
  }
  fit[fit>1]<-1
  fit[fit<0]<-0
  if(type=="probs"|type=="prob"){
    return(cbind(1-fit,fit))
  }
  return(list(class=ty,probs=cbind(1-fit,fit)))
}
