"predict.svmpath" <-
  function(object,newx,lambda,type=c("function","class","alpha","margin"),...){
    type<-match.arg(type)
    oalpha<-cbind(object$alpha)# in case its not a matrix
    oalpha0<-object$alpha0
    if(missing(lambda)){
      lambda<-object$lambda
      alpha<-oalpha
      alpha0<-oalpha0
    }
    else if (length(object$lambda)==1){
    warning("This model has only one lambda value; predictions at it are returned")
    nlam=length(lambda)
    lambda=rep(object$lambda,nlam)
    alpha=oalpha[rep(1,nlam),,drop=FALSE]
    alpha0=rep(oalpha0,nlam)
  }
  else{
      olambda<-object$lambda
      nalpha<-length(object$y)
      minl<-min(olambda);maxl<-max(olambda)
      anysmaller<-seq(along=lambda)[lambda<minl]
      if(length(anysmaller))lambda[anysmaller]<-minl
      lmax<-max(lambda)
      if(lmax>maxl){# we need to modify our alphas
        alpha00<-object$alpha00
        maxl<-lmax
        oalpha<-cbind(oalpha[,1],oalpha)
        oalpha0<-c(alpha00["slope"]*lmax+alpha00["intercept"],oalpha0)
        olambda<-c(maxl,olambda)
      }
      lfrac<-(lambda-minl)/(maxl-minl)
      olfrac<-(olambda-minl)/(maxl-minl)
      coord<-approx(olfrac,seq(olambda),lfrac)$y
      left<-floor(coord);right<-ceiling(coord)
      alpha <- outer(rep(1,nalpha),olfrac[right] - lfrac) * oalpha[,left , drop = FALSE] + 
                outer(rep(1,nalpha),lfrac - olfrac[left]) * oalpha[,right , drop = FALSE]
      alpha<-scale(alpha,  FALSE, olfrac[right] - olfrac[left])
      alpha[,left == right] <- oalpha[,left[left == right] ]
      alpha0 <- ((olfrac[right] - lfrac) * oalpha0[left] + 
                 (lfrac - olfrac[left]) * oalpha0[right])/(olfrac[right] - 
                                                           olfrac[left])
      alpha0[left == right] <- oalpha0[left[left == right] ]
      }
    if(type=="alpha"){
      attr(alpha,"scaled:scale")<-NULL
      object<-list(alpha0=alpha0,alpha=drop(alpha),lambda=lambda)
    }
    else{
      if(missing(newx))newx<-object$x
      K<-object$kernel(newx,object$x,object$param.kernel)
      fit<-K%*%(alpha*object$y)
      if(type=="margin"){
        margin=(alpha*object$y)*fit
        margin=apply(margin,2,sum)
        margin=lambda/sqrt(margin)
      }
      fit<- scale(fit,-alpha0,lambda)
      attr(fit,"scaled:center")<-NULL
      attr(fit,"scaled:scale")<-NULL
    }
    switch(type,
           "function"=fit,
           "class"=sign(fit),
           alpha=object,
           margin=margin
           )
  }
