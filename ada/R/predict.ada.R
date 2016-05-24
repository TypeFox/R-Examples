"predict.ada" <-
function(object,newdata=NULL,type=c("vector","probs","both","F"),n.iter=NULL,...){
  if(!inherits(object,"ada")){
    stop("Error:  Object is not of class ada")
  }
  if(missing(type)){
    type="vector"
  }
  if(type=="probs")
    type="prob"
  if(type!="vector" & type!="prob" & type!="both" & type!="F"){
    warning(paste("type=",type," is undefined:  default is 'vector'..  "))
    type="vector"
  }
  if(missing(newdata)){
    stop("Error: Arguement newdata is missing and must be supplied\n")
  }
  iter=object$iter
  if(!is.null(n.iter)){
    if(n.iter<iter)
      iter=n.iter
  }
  lev<-levels(as.factor(object$actual))
  const<-2
 
  f<-object$model$lossObj$predict.type
  tmp=sapply(1:iter,function(i)f(f=object$model$trees[[i]],dat=newdata))
  tmp=t(t(tmp)*object$model$alpha[1:iter])
  tmp<-apply(tmp,1,sum)
  fit<-as.vector(sign(tmp))
##fit<-as.factor(fit)
  a1=(fit==0)
  if(sum(a1)>0)
    fit[a1]<-sample(c(-1,1),sum(a1),TRUE,c(.5,.5))
  ind<-as.numeric(factor(fit,levels=c(-1,1)))
  fit<-factor(lev[ind],levels=lev)
	
##attr(fit,"levels")<-lev
  if(type=="vector")
    return(fit)
   cal<-function(x,const){
     if(x>0)
       return(c(exp(-const*x),1))
     return(c(1,exp(const*x)))
   }
  probs<-t(sapply(tmp,cal,const=const))
  probs<-probs/apply(probs,1,sum)
  if(type=="prob")
    return(probs)
  if(type=="F")
    return(tmp)
  return(list(class=fit,probs=probs))
}

