"update.ada" <-
function(object,x,y,test.x,test.y=NULL, n.iter,...){
  if(class(object)!="ada"){
    stop("Update routine is for an ada object only")
  }
  if(missing(y) | missing(x)){
    stop("This procedure requires a response and a set of variables")
  }
  if(missing(n.iter)){
	stop("The new number of iterations must be specified")
  } 
  if(n.iter<=object$iter){
    return(object)
  }
  lev=levels(as.factor(y))
  x=as.data.frame(x)
  test.dat=FALSE
  if(!is.null(test.y)){
    tlev=NULL
    if(!missing(test.x)){
      tlev=levels(as.factor(test.y))
      test.x=as.data.frame(test.x)
      test.dat=TRUE
    }
    if(length(tlev)<1)
      warning("test.x must be the testing data and the response must have 2 levels")
  }else{
    test.x=NULL
  }
  ny<-c(-1,1)[as.numeric(as.factor(y))]
  if(test.dat)
    test.y<-c(-1,1)[as.numeric(as.factor(test.y))]
  control=as.list(object$model$trees[[1]]$control)
  na.action=object$na.action
  lossObj=object$model$lossObj
  result =ada.machine(x,ny,test.x,test.y,n.iter,lossObj=lossObj,oldObj=object,control=control,na.action=na.action)
  g=as.factor(lev[as.numeric(as.factor(sign(result$F[[1]])))])
  tab=table(as.factor(y),g,dnn=c("True value","Final Prediction"))
  obj=structure(list(model=result,fit=g,call=object$call,confusion=tab,iter=n.iter,
    actual=as.vector(y),nu=object$nu,dim=object$dim,bag.frac=object$bag.frac,names=object$names),class="ada")
  obj
}

