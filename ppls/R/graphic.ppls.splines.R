graphic.ppls.splines=function(X,y,lambda=NULL,add.data=FALSE,select=FALSE,ncomp=1,deg=3,order=2,nknot=NULL,reduce.knots=FALSE,kernel=TRUE,window.size=c(3,3)){
    p<-ncol(X)
    ntest<-300 # number of test examples
    Xtest<-matrix(,ntest,p)
    for (i in 1:p){
        Xtest[,i]=seq(min(X[,i]),max(X[,i]),length=300) 
    }
    # transform training and test data
    Z<-X2s(X,Xtest,deg=deg,nknot=nknot,reduce.knots=reduce.knots)
    Ztrain<-Z$Z
    Ztest<-Z$Ztest
    sizeZ<-Z$sizeZ
    P <- lambda*Penalty.matrix(m =sizeZ,order=order)
    blocks=c()
    for (b in 1:length(sizeZ)) {
            blocks = c(blocks, rep(b, sizeZ[b]))
        }

    ppls.object<-penalized.pls(Ztrain,y,P=P,ncomp=ncomp,select=select,kernel=kernel,blocks=blocks)
    ppls.coefficients<-ppls.object$coefficients[,ncomp]
    Ytest<-matrix(,ntest,ncol(X)) # prediction for each additive component
    for (i in 1:ncol(X)){
        start<-cumsum(c(0,sizeZ))[i]+1 # start of the ith block
        end<-cumsum(sizeZ)[i] # end of the ith block
        Ytest[,i]<-Ztest[,start:end]%*%ppls.coefficients[start:end]
    }
    # plot the predicted functions
    par(mfrow=window.size)
    for (i in 1:p){
        plot(Xtest[,i],Ytest[,i],type="l",lwd=3,xlab="x",ylab="y",main=i,col="blue")
        if (add.data==TRUE){
        lines(X[,i],scale(y,scale=FALSE),type="p",lwd=2)
        lines(Xtest[,i],Ytest[,i],type="l",lwd=3,col="blue")
        }
    }  
    return(ppls.coefficients)  

}