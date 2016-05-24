gelman.prior<-function(formula, data, scale=1, intercept=scale, singular.ok=FALSE){
  X1<-model.matrix(formula, data)
  if(singular.ok==FALSE){
     sing.rm<-lm(rnorm(nrow(X1))~X1-1)
     sing.rm<-which(is.na(sing.rm$coef))
     if(length(sing.rm)>0){
       warning("some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!")
     }
   }	

  X2<-get_all_vars(formula, data)
  X2<-as.data.frame(lapply(X2, function(x){if(is.numeric(x)){scale(x, scale=sd(x)*2*(length(x)-1)/length(x))}else{x}}))
  X2<-model.matrix(formula, data=X2)
  if(all(X2[,1]==1)){
    X2[,-1]<-apply(X2[,-1,drop=FALSE], 2, function(x){if(any(!x%in%c(0,1))){x}else{scale(x, center=sum(x)/length(x), scale=1)}})
  }else{
    X2<-apply(X2, 2, function(x){if(any(!x%in%c(0,1))){x}else{scale(x, center=sum(x)/length(x), scale=1)}})
  }
  if(length(sing.rm)>0){
    X1<-X1[,-sing.rm]
    X2<-X2[,-sing.rm]
  }
  P<-solve(t(X1)%*%X1, t(X1)%*%X2)
  I<-diag(nrow(P))*scale^2
  I[1,1]<-intercept^2
  P%*%I%*%t(P)
}

