"update.maboost" <-
  function(object,x,y,test.x,test.y=NULL, n.iter,...){
    if(class(object)!="maboost"){
      stop("Update routine is for an maboost object only")
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
    K=length(lev)
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
    if(K==2){
      multiclass=FALSE;
    ny<-c(-1,1)[as.numeric(as.factor(y))]
    if(test.dat)
      test.y<-c(-1,1)[as.numeric(as.factor(test.y))]
    }else{
      multiclass=TRUE;
      ny<-factor(y,levels=1:K)
      if(test.dat)
        test.y<-factor(test.y,levels=1:K)
      
    }
    control=as.list(object$model$trees[[1]]$control)
    na.action=object$na.action
    lossObj=object$model$lossObj
    maxmargin=object$maxmargin
    smooth=object$smooth
    smoothfactor=object$smoothfactor
    nu=object$nu
    bag.frac=object$bag.frac
    random.feature=object$random.feature
    random.cost=object$random.cost
    Ctree=object$Ctree
    cl=object$call
    if(multiclass){
      result =maboost.machine.multi(x,ny,test.x,test.y,n.iter,maxmargin,smooth,smoothfactor,nu,bag.frac,random.feature,random.cost,lossObj,oldObj=NULL,na.action=na.action,...)
      
      g = factor ( apply(result$F[[1]],1,which.max),levels=1:K )
      levels(g)=lev;
      
    }else{
      
      result =maboost.machine.bin(x,ny,test.x,test.y,n.iter,maxmargin,smooth,smoothfactor,nu,bag.frac,random.feature,random.cost,lossObj,oldObj=NULL,na.action=na.action,...)
      g  =  factor(sign(result$F[[1]]),levels=c(-1,1));
      levels(g)=lev;
      
      
      
    }
    
    
    tab=table(as.factor(y),g)
    nm<-1:(dim(x)[2])
    if(is.data.frame(x)){
      nm=names(x)
    } 
    obj=structure(list(model=result,fit=g,call=cl,confusion=tab,iter=n.iter,
                       actual=as.factor(y),nu=nu,dim=dim(x),names=nm,bag.frac=bag.frac
                       ,na.action=na.action,
                       Ctree=Ctree,maxmargin=maxmargin
                       ,smooth=smooth,
                       smoothfactor=smoothfactor,
                       random.feature=random.feature,
                       random.cost=random.cost)
                  ,class="maboost")
    obj
  }
