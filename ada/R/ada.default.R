"ada.default" <-
function(x,y,test.x,test.y=NULL,loss=c("exponential","logistic"),
         type=c("discrete","real","gentle"),iter=50, nu=0.1, bag.frac=0.5,
         model.coef=TRUE,bag.shift=FALSE,max.iter=20,delta=10^(-10),verbose=FALSE,
         ...,na.action=na.rpart){
  cl<-match.call(expand.dots=TRUE)
  cl[[1]]<-as.name("ada")
  
  type = match.arg(type)
  if(missing(type)){
    type="discrete"
  }
  if(missing(loss))
    loss="exponential"
  loss=tolower(loss)
  if(loss=="e"|loss=="ada"){
    loss="exponential"
  }
  if(loss=="l"|loss=="l2"){
    loss="logistic"
  }
  if(missing(y) | missing(x)){
    stop("This procedure requires a response and a set of variables")
  }
  lev=levels(as.factor(y))
  if(length(lev)!=2)
    stop("Currently this procedure can not directly handle > 2 class response\n")
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
  if(length(lev)>2){
    stop(paste("Error:  At this time ",cl[[3]]," must have 2 levels",sep=""))
  }
  extraArgs=list(...)
  if (length(extraArgs)) {
    arg <- names(formals(rpart))
    indx <- match(names(extraArgs), arg, nomatch = 0)
    if (any(indx == 0)) 
      stop(paste("Error:  Argument", names(extraArgs)[indx == 0], 
                 "not matched"))
  }
  ny<-c(-1,1)[as.numeric(as.factor(y))]
  if(test.dat)
    test.y<-c(-1,1)[as.numeric(as.factor(test.y))]

  ### Set Up Predictions for each boosting type
  method="class"
  if(type=="discrete"){
	predict.type<-function(f,dat){ 
	 		p <- predict(f,newdata=dat,type="class") 
	 		a=as.numeric(levels(p))[p] 
	 		a 
 	}
  }
  if(type=="real"){
    predict.type<-function(f,dat){
      a=predict(f,newdata=dat,type="prob")[,2]
      ind<-a==1|(1-a)==1
      if(length(ind)>1){
        a[ind]=(1-a[ind])*0.0001+a[ind]*.9999
      }
      f=0.5*log(a/(1-a))
      if(sum(is.nan(f))>0)
        f[is.nan(f)]=0.51
      f
    }
  }
  if(type=="gentle"){
    predict.type<-function(f,dat){
      predict(f,newdata=dat)
    }
    method="anova"
  }

  ### Set up coeficients
  wfun<-function(yf,f){
    exp(-yf*f)
  }
  if(loss=="logistic"){
    wfun<-function(yf,f){
      a=exp(-yf*f)
      a/(1+a)
    }
  }
  coefs=function(wts,eta,yf,alp){
    1
  }
  if(model.coef){
    if(loss=="exponential"){
      if(type=="discrete"){ ### eta is assumed = err
        coefs=function(wts,eta,yf,alp){
          alp
        }
      }else{
        coefs=function(wts,eta,yf,alp){
          alpst=alp
          for(i in 1:max.iter){
            alp0=alp
            pval=wts*exp(-yf*eta*alp)
            m11<-sum(eta*pval*eta)
            m12<-sum(yf*pval*eta)
            if(m11==0){
              alp=alpst
              break
            }
            alp=alp+m12/m11
            a1=(alp-alp0)^2
            if(a1<delta){
              break
            }
          }
          if(verbose)
            cat("FINAL: iter=",i," rate=",a1,"\n")
          alp
        }
      }
    }else{
      coefs=function(wts,eta,yf,alp){
        alpst=alp
        for(i in 1:max.iter){
          alp0=alp
          pval=wts*exp(-yf*eta*alp)
          pval=1-1/(1+pval)
          m11<-sum(eta*pval*eta)
          m12<-sum(yf*pval*eta)
          if(m11==0){
            alp=alpst
            break
          }
          alp=alp+m12/m11
          a1=(alp-alp0)^2
          if(a1<delta){
            break
          }
          alp=alp+m12/m11
        }
        if(verbose)
          cat("FINAL: iter=",i," rate=",a1,"\n")
        2*alp
      }
    }
  }
  lossObj=list(coefs=coefs,wfun=wfun,predict.type=predict.type,method=method,type=type,loss=loss,shift=bag.shift)
  result =ada.machine(x,ny,test.x,test.y,iter,nu,bag.frac,lossObj,oldObj=NULL,na.action=na.action,...)
  g=as.factor(lev[as.numeric(as.factor(sign(result$F[[1]])))])
  tab=table(as.factor(y),g,dnn=c("True value","Final Prediction"))
  nm<-1:(dim(x)[2])
  if(is.data.frame(x)){
    nm=names(x)
  } 
  obj=structure(list(model=result,fit=g,call=cl,confusion=tab,iter=iter,
    actual=as.vector(y),nu=nu,dim=dim(x),names=nm,bag.frac=bag.frac,na.action=na.action),class="ada")
  obj
}

