GORMC <-
function(survfun=formula(data),curefun=formula(data),data=parent.frame(),
                 r=0,n.int=5,order=3,max.iter=1000,cov.rate=.001){
  sdata<-data
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  temp <- c("", "formula", "data", "na.action")
  mf <- mf[match(temp, names(mf), nomatch = 0)]
  mf[[1]] <- as.name("model.frame")
  avars<-all.vars(survfun)
  survfun1<-as.formula(paste("Surv(",avars[1],",",avars[2],",type='interval2')~",paste(avars[-c(1,2)],collapse="+"),sep=""))

  temp.x<-terms(survfun1, data = sdata)
  temp.z<-terms(curefun, data = sdata)
  mf$formula<-temp.x
  mf <- eval(mf, envir=parent.frame())
  Y <- model.extract(mf, "response")
  attr(temp.x, "intercept") <- 0
  attr(temp.z, "intercept") <- 1
  Xp <- model.matrix(temp.x,sdata)
  colnames(Xp)<-paste("x",1:ncol(Xp),sep="")
  Zp <- model.matrix(temp.z,sdata)
  colnames(Zp)<-paste("z",c(1:ncol(Zp))-1,sep="")

  sdata$d1<-as.numeric(Y[,3]==2)
  sdata$d2<-as.numeric(Y[,3]==3)
  sdata$d3<-as.numeric(Y[,3]==0) 
  sdata$Li<-Y[,1]
  sdata$Li[sdata$d1==1]<-0
  sdata$Ri<-Y[,2]
  sdata$Ri[sdata$d1==1]<-Y[sdata$d1==1,1]
  sdata$Ri[sdata$d3==1]<-NA
   
   L<-n.int+order
   P<-ncol(Xp)
   M<-ncol(Zp)
   N<-nrow(sdata)

   # Estimates
   temp.est<-try(EM.Iter(sdata=sdata,Xp=Xp,Zp=Zp,r=r,n.int=n.int,order=order,max.iter=max.iter,cov.rate=cov.rate),silent=TRUE)
   if(is.character(temp.est)) cat("ERROR: the EM algorithm does not converge.\n")
   if(!is.character(temp.est)){
   # Variance
   temp.hess<-try(VCOV.calc(object=temp.est),silent=TRUE)
     if(is.character(temp.est)) cat("ERROR: the variance-covariance matrix is not available.\n")
     if(!is.character(temp.est) & temp.hess$hess==1) cat("NOTE: the QR decomposition is applied to solve the hessian matrix.\n")
     if(!is.character(temp.est) & temp.hess$hess==2) cat("NOTE: the g-inverce is applied to the hessian matrix.\n")
     if(!is.character(temp.est) & temp.hess$hess==3) cat("NOTE: the hessian matrix is obtained from numerical methods.\n")
   }
output<-list(survfun=survfun1,curefun=curefun,ParEst=temp.est,ParVcov=temp.hess,call=match.call())
class(output)<-"GORMC"
output$mdata<-sdata
return(output)
}
