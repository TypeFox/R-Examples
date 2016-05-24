"addtest" <-
function(x,test.x,test.y,...){
  if(!inherits(x,"ada")){
    stop("Error:  Object is not of class ada")
  }
  if(missing(test.y)){
    stop("This funciton needs a tesing response")
  }
  if(missing(test.x)){
    stop("This function needs testing data")
  }
  kapstat<-function(tab=diag(2) ){
    if(dim(tab)[1]==1){
      return(0)
    }
    if(dim(tab)[1]==dim(tab)[2]){
      rs<-apply(tab,2,sum)
      cs<-apply(tab,1,sum)
      N<-sum(rs)
      E<-sum(rs*cs)/N^2
      O<-sum(diag(tab))/N
      return( (O-E)/(1-E) )
    }else{
      return(NA)
    }
  }
  iter=x$iter
  lev=levels(as.factor(test.y))
  if(length(lev)>2){
    stop("Error The response must have 2 levels")
  }
  y=c(-1,1)[as.numeric(as.factor(test.y))]
  nt<-"vector"
  test.errs<-test.kaps<-rep(0,iter)
  test.x=as.data.frame(test.x)

  f<-x$model$lossObj$predict.type
  tmp=sapply(1:iter,function(i)f(f=x$model$trees[[i]],dat=test.x))
  tmp=t(t(tmp)*x$model$alpha)
  Fm=rep(0,length(y))
  for(m in 1:iter){
    Fm<-Fm+tmp[,m]
    tab<-matrix(table(sign(Fm),y),nrow=2,ncol=2)
    test.errs[m]<-1-sum(diag(tab))/sum(tab)
    test.kaps[m]<-1-kapstat(tab)
  }
  x$model$errs=cbind(x$model$errs,test.errs,test.kaps)
  n1<-length(x$model$F)
  x$model$F[[n1+1]]<-Fm
  x
}

