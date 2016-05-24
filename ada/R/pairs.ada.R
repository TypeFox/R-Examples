"pairs.ada" <-
function(x,train.data=NULL,vars=NULL,maxvar=10,test.x=NULL,test.y=NULL,test.only=FALSE,col=c(2,4),pch=c(1,2),...){
  if(class(x)!="ada"){
    stop("Error:  Object is not of class ada")
  }
  if(is.null(train.data) & is.null(test.x)){
    stop("Note:  The train.data must correspond to the predictors used in adaboost")
  }
  if(maxvar<2){
    warning("Note: maxvar must be greater than 2.  So this variable will be ignored...")
    maxvar=10
  }
  iter<-x$iter
  y<-as.numeric(as.factor(as.vector(x$actual)))+1
  lev<-levels(as.factor(y))
  if(!is.null(train.data))
    ptrain<-predict(x,newdata=train.data,type="both")
  ptest<-list(class=NULL,probs=NULL)
  if(!is.null(test.x)){
    if(!is.null(train.data)){
      if(dim(train.data)[2]!= dim(test.x)[2]){
        stop("Error:  The test data set must have the same number of variables as the train.data")
      }
    }
    ptest<-predict(x,test.x,type="both")
    if(is.null(test.y)){
      test.y<-rep(1,dim(test.x)[1])
      y<-c(y,as.numeric(test.y))
    }else{
      test.y=as.numeric(lev[as.numeric(as.factor(test.y))])
    }
  }
  ptrain$class<-c(ptrain$class,ptest$class)
  ptrain$probs<-rbind(ptrain$probs,ptest$probs)
  if(test.only & !is.null(test.x)){
    ptrain$class<-c(ptest$class)
    ptrain$probs<-rbind(ptest$probs)
    y<-as.numeric(test.y)
    train.data=NULL
  }
  tr<-y
  mat<-ptrain
  train.data<-as.data.frame(train.data)
  if(is.null(vars)&(!test.only)){
    if(dim(train.data)[2]>maxvar){
      if(!is.data.frame(train.data)){
        stop(" train.data must be of type data.frame")
      }
      vars<-names(varplot(x,plot.it=FALSE,type="scores",max.var.show=maxvar))
    }else{
      vars<-1:(x$dim[2])
    }
  }
  if(test.only & is.null(vars)){
    if(dim(test.x)[2]>maxvar){
      if(!is.data.frame(test.x)){
        stop(" test.x must be of type data.frame")
      }
      vars<-names(varplot(x,plot.it=FALSE,type="scores",max.var.show=maxvar))
      ##match(vars,val$names)
    }else{
      vars<-1:(x$dim[2])
    }
  }
 panel.up<-function(x,y,...){
   val<-as.numeric(as.factor(tr))
   points(x,y,col=col[val],pch=pch[val])
 }
 panel.low<-function(x,y,...){
   val<-as.numeric(as.factor(mat$class))
   points(x,y,col=col[val],cex=apply(mat$probs,1,max),pch=pch[val])
 }
 pairs(as.matrix(rbind(train.data,test.x))[,vars],lower.panel=panel.low,upper.panel=panel.up,...)
}

