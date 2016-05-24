ridge.cv<-
function(X,y,lambda=NULL,scale=TRUE,k=10,plot.it=FALSE){
        if (is.vector(X)==TRUE){
            X<-matrix(X,ncol=1)
        }
    if (is.null(lambda)==TRUE){
        ss<-seq(-10,-1,length=1000)
        ss<-10^ss
        n<-nrow(X)
        nn<-n- floor(n/k)
        lambda<-ss*nn*ncol(X)
    }
    cv<-rep(0,length(lambda))
    n<-nrow(X)
    all.folds <- split(sample(1:n), rep(1:k,length=n))
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain=X[-omit,,drop=FALSE]
        ytrain=y[-omit]
        Xtest=X[omit,,drop=FALSE]
        ytest=y[omit]
        if ((is.vector(X)==TRUE)| (ncol(X)==1)){
              xtrain<-as.vector(Xtrain)
              coef.ll<-lm.ridge.univariate(xtrain,ytrain,lambda=lambda,scale=scale)
        }
        else{
        ll<-lm.ridge(ytrain~Xtrain,scale=scale,lambda=lambda)
        coef.ll<-coef(ll)
    }
        res<-matrix(,length(ytest),length(lambda))
        pred<-t(matrix(coef.ll[,1],nrow=length(lambda),ncol=length(ytest))) + Xtest%*%t(coef.ll[,-1])
        res<-pred-matrix(ytest,nrow=length(ytest),ncol=length(lambda))
        cv<-cv+apply(res^2,2,sum)
        
    }
    cv<-cv/n
    lambda.opt<-lambda[which.min(cv)]
    if (plot.it==TRUE){
        plot(lambda,cv,type="l")
    }
    if ((is.vector(X)==TRUE)| (ncol(X)==1)){
      x<-as.vector(X)
     coefficients<- as.vector(lm.ridge.univariate(x,y,scale=scale,lambda=lambda.opt))
      
    }
    else{
    rr<-lm.ridge(y~X,scale=scale,lambda=lambda.opt)
    coefficients<-coef(rr)
    }
    intercept<-coefficients[1]
    coefficients<-coefficients[-1]
    return(list(intercept=intercept,coefficients=coefficients,lambda.opt=lambda.opt))
}
