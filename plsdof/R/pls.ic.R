pls.ic=function (X, y, m = min(ncol(X),nrow(X)-1),criterion="bic",naive=FALSE,use.kernel=FALSE,compute.jacobian=FALSE,verbose=TRUE) {
    m.crash=NA
    n <- nrow(X)
    DoF.max = min(n-1, ncol(X)+1)
    compute.DoF=TRUE
    if (naive==TRUE){
        compute.DoF=FALSE
        compute.jacobian=FALSE
    }
    pls.object <- pls.model(X,y,m,compute.DoF=compute.DoF,compute.jacobian=compute.jacobian,use.kernel=use.kernel)
    RSS <- pls.object$RSS
    yhat <- pls.object$yhat
    sigmahat <- (pls.object$sigmahat)
        DoF <- pls.object$DoF
    if (min(DoF)<=0){
        sign.DoF<-sign(DoF)
        sign.DoF[sign.DoF==0]=-1
        dummy<-(1:(m+1))[sign.DoF==-1]
        mini<-min(dummy)-1
        m<-mini-1
        m.crash<-m+1
        if (verbose==TRUE){
        if (compute.jacobian==TRUE){
        cat(paste("Negative DoF for jacobian. Setting maximal number of components to ",m,".\n"))
        }
        if (compute.jacobian==FALSE){
        cat(paste("Negative DoF. Setting maximal number of components to ",m,".\n"))
        }
        }
        #cat(paste("DoF for ",m," components: ",DoF[1,m+1],"\n"))
        DoF[(m:length(DoF))]=Inf
    }
    if (min(DoF)>0){
    ic <- information.criteria(RSS, DoF, yhat = yhat, sigmahat = sigmahat, 
        n, criterion=criterion)
    #DoF = ic$DoF
    m.opt <- ic$par-1
        score<-ic$score
}
coefficients<-pls.object$coefficients[,m.opt]
intercept<-pls.object$intercept[m.opt]
covariance<-pls.object$covariance
if (compute.jacobian==TRUE){
    covariance=covariance[m.opt+1,,]
}
outlist=list(DoF = DoF, m.opt = m.opt,sigmahat=sigmahat,m.crash=m.crash,intercept=intercept,coefficients=coefficients,covariance=covariance)
class(outlist)="plsdof"    
    return(outlist)
}
