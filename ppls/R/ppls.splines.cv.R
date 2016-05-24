`ppls.splines.cv` <-
function(X,y,lambda=1,ncomp=NULL,degree=3,order=2,nknot=NULL,k=5,kernel=FALSE,scale=FALSE,reduce.knots=FALSE,select=FALSE){

  n<-nrow(X)

  p<-ncol(X)

  if (is.null(ncomp)) ncomp=min(n-1,p)
  
  lambda=as.vector(lambda)

 
  all.folds <- split(sample(1:n), rep(1:k,length=n))
  
  # ensure that ncomp does not exceed the sample size on the cv splits
  ntrain=vector(length=k)
  for (i in 1:k){
    ntrain[i]=n-length(all.folds[[i]])
  }
  ntrain.min=min(ntrain)
  ncomp=min(ncomp,ntrain.min-1)
  #
   error.cv=matrix(0,length(lambda),ncomp)

  for (i in seq(k)){
  
        omit <- all.folds[[i]]
        Xtrain=X[-omit,,drop=FALSE]
        ytrain=y[-omit]
        Xtest=X[omit,,drop=FALSE]
        ytest=y[omit]
                
        Z<-X2s(Xtrain,Xtest,degree,nknot,reduce.knots=reduce.knots)

        Ztrain=Z$Z

        Ztest<-Z$Ztest

        P<-Penalty.matrix(m=Z$sizeZ,order=order)
        blocks=c()
        for (b in 1:length(Z$sizeZ)){
            blocks=c(blocks,rep(b,Z$sizeZ[b]))
        }
    
        for (j in 1:length(lambda)){
            penpls=penalized.pls(Ztrain,ytrain,lambda[j]*P,ncomp,kernel,blocks=blocks,select=select,scale=scale)
            error.cv[j,]=error.cv[j,]+ length(ytest)*(new.penalized.pls(penpls,Ztest,ytest)$mse)
        }
  
  }
  #cat(paste("cv completed \n"))
  error.cv=error.cv/n
  value1=apply(error.cv,1,min)
  lambda.opt=lambda[which.min(value1)]
  ncomp.opt=which.min(error.cv[lambda==lambda.opt,])
  min.ppls=min(value1)
  return(list(error.cv=error.cv,min.ppls=min.ppls,lambda.opt=lambda.opt,ncomp.opt=ncomp.opt))
}
