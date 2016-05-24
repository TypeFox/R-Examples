per.variance<- function(object){  
  per <- NULL  
  if((prod(object$keepY==rep(dim(object$Y)[2],object$ncomp))==1)|is.null(object$keepY)) {
  for (h in 1:object$ncomp)  {
    per <- c(per,mean(apply(matrix(object$variates$X[,h],ncol=1)%*%matrix(t(object$mat.d)[h,],nrow=1),MARGIN=2,FUN=var)/apply(object$Y,MARGIN=2,FUN=var)))
  }
  }else{
    for (h in 1:object$ncomp)  {
      if(class(object)[1]=="sPLS") {
      YY <- object$Y[,select.spls(object)$select.Y.total]
      ind <- select.spls(object)$select.Y.total} else{
      YY <- object$Y[,select.sgpls(object)$select.Y.total]
      ind <- select.sgpls(object)$select.Y.total
      }
      per <- c(per,mean(apply(matrix(object$variates$X[,h],ncol=1)%*%matrix(t(object$mat.d)[h,ind],nrow=1),MARGIN=2,FUN=var)/apply(YY,MARGIN=2,FUN=var)))
    } 
    
    
  }
  result <- list(perX=per,cum.perX=cumsum(per))
  return(result)
}