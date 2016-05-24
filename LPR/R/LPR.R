#LPR
LPR<-function(x, y, lambda2, B, type.boot=c("residual","paired"), alpha = 0.05){
   mLS<-TRUE
   type.lambda<-"cv"
   x<- as.matrix(x)
   y<- as.numeric(y)
   n<- dim(x)[1]
   p<- dim(x)[2]
  
  if(missing(B)){
    B<-500
  }
  if(missing(lambda2)){
    lambda2<-1/n
  }
  
 
  selectset<-rep(0,p)
  Beta<-rep(0,p)
  Beta.LPR<-rep(0,p) 
  
  
    globalfit<-glmnet(x, y, standardize=FALSE,intercept=FALSE)	
    lambda<-globalfit$lambda
    cvfit<-escv.glmnet(x, y, lambda=lambda, nfolds=5, tau=0, mLS=mLS, standardize=FALSE, 
                       intercept=FALSE)
    
    if(type.lambda=="cv"){
      lambda.opt<-cvfit$lambda.cv
    }
    if(type.lambda=="escv"){
      lambda.opt<-cvfit$lambda.escv
    }    
    
    fitlasso <- predict(globalfit, type = "coefficients", s = lambda.opt)
    fit.value<- predict(globalfit, newx=x, s = lambda.opt)
    betalasso<-fitlasso[-1]
    Beta<-betalasso
    selectvar<-betalasso!=0
    selectset<-as.numeric(selectvar)
   
    PR.obj<-PartRidge(x=x, y=y, lambda2=lambda2, varset=selectvar, standardize=FALSE, intercept=FALSE)
    Beta.LPR<-PR.obj$beta			
    
    if(type.boot=="residual"){
      fit<-x%*%Beta
      residual<-y-fit
      residual_center<-residual-mean(residual)
      Beta.boot.LPR<-matrix(0,nrow=B,ncol=p)
      
      for(i in 1:B){
        resam<-sample(1:n,n,replace=TRUE)   			
        ystar<-fit+residual_center[resam]
        if(!mLS){
          boot.obj<-mylasso(x, ystar, lambda=lambda.opt, standardize=FALSE, intercept=FALSE)
        } 
        if(mLS){
          boot.obj<-mylassomls(x, ystar, lambda=lambda.opt, standardize=FALSE, intercept=FALSE)
        }	
        boot.selectvar<-boot.obj$beta!=0			
        boot.PR.obj<-PartRidge(x=x, y=ystar, lambda2=lambda2, varset=boot.selectvar, standardize=FALSE, intercept=FALSE)
        Beta.boot.LPR[i,]<-boot.PR.obj$beta
      }
      
      interval.LPR<-ci(Beta.LPR, Beta.boot.LPR, alpha=alpha, type="basic")
      }
    
     if(type.boot=="paired"){
      Beta.boot<-matrix(0,nrow=B,ncol=p)
      Beta.boot.LPR<-matrix(0,nrow=B,ncol=p)
      
      for(i in 1:B){
        resam<-sample(1:n,n,replace=TRUE)			
        rx<-x[resam,]
        ry<-y[resam]				
        if(!mLS){
          boot.obj<-mylasso(rx, ry, lambda=lambda.opt, standardize=FALSE, intercept=FALSE)
        } 
        if(mLS){
          boot.obj<-mylassomls(rx, ry, lambda=lambda.opt, standardize=FALSE, intercept=FALSE)
        }	
        Beta.boot[i,]<-boot.obj$beta
        boot.selectvar<-boot.obj$beta!=0			
        boot.PR.obj<-PartRidge(x=rx, y=ry, lambda2=lambda2, varset=boot.selectvar, standardize=FALSE, intercept=FALSE)
        Beta.boot.LPR[i,]<-boot.PR.obj$beta
      }
      
    
      interval.LPR<-ci(Beta.LPR, Beta.boot.LPR, alpha=alpha, type="quantile")
       }			
    
  object<-list(lambda.opt=lambda.opt, Beta=Beta, Beta.LPR=Beta.LPR,  
               n=n, p=p, selectset=selectset, interval.LPR=interval.LPR
               )	 
  object
}






