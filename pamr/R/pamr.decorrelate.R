pamr.decorrelate<-  function (x, adjusting.predictors, xtest=NULL, adjusting.predictors.test=NULL){

foo<- lm(t(x)~., adjusting.predictors)
x.adj=t(foo$res)
xtest.adj=NULL

if(!is.null(adjusting.predictors.test)){
   if(is.null(xtest)){
   stop("xtest must be non-null if adjusting.predictors.test is non-null")
  }
    temp=t(predict(foo,adjusting.predictors.test))
    xtest.adj=xtest-temp
}
return(list(x.adj=x.adj,xtest.adj=xtest.adj))
}
