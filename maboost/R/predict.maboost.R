"predict.maboost" <-
  function(object,newdata=NULL,type=c("class","prob","both","F"),n.iter=NULL,...){
    if(!inherits(object,"maboost")){
      stop("Error:  Object is not of class maboost")
    }
    if(missing(type)){
      type="class"
    }
    if(type=="probs")
      type="prob"
    if(type!="class" & type!="prob" & type!="both" & type!="F"){
      warning(paste("type=",type," is undefined:  default is 'class'..  "))
      type="class"
    }
    if(missing(newdata)){
      stop("Error: Arguement newdata is missing and must be supplied\n")
    }
    iter=object$iter
    if(!is.null(n.iter)){
      if(n.iter<iter)
        iter=n.iter
    }
    repmat = function(X,m,n){
      ##R equivalent of repmat (matlab)
      mx = length(X)
      nx = 1
      matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
    }
    lev<-levels(as.factor(object$actual))
    K=length(lev);
    if (length(lev)>2){
      M=dim(newdata)[1]
      
      pred<-object$model$lossObj$predict.type
      tmp = matrix(0,M,length(lev))
      for ( i in 1 : iter){
        ff = pred(f=object$model$trees[[i]],dat=newdata)
        tmp[cbind( (1:M),ff )] = tmp[cbind( (1:M),ff )] + object$model$alpha[i];
      }
      fit = apply(tmp,1,which.max)
      fit<-factor(fit,levels=1:K)
     
      
      ##fit<-as.factor(fit)
      levels(fit)=lev
         
      ##attr(fit,"levels")<-lev
      if(type=="class")
        return(fit)
      
      probs<-tmp/repmat(apply(tmp,1,sum),1,K)
      colnames(probs)=lev;
      if(type=="prob"){
        
        return(probs)
      }
      if(type=="F"){
        return(tmp)
      }
      
    }else{
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
      if(type=="class")
        return(fit)
      cal<-function(x,const){
        if(x>0)
          return(c(exp(-const*x),1))
        return(c(1,exp(const*x)))
      }
      probs<-t(sapply(tmp,cal,const=const))
      probs<-probs/apply(probs,1,sum)
      colnames(probs)=lev;
      if(type=="prob")
        return(probs)
      if(type=="F")
        return(tmp)
      
      
    }
    return(list(class=fit,probs=probs))
    
  }
