ridge.cv<-
function(X,y,lambda=NULL,scale=TRUE,k=10,plot.it=FALSE,groups=NULL,method.cor="pearson",compute.jackknife=TRUE){
    if (is.null(lambda)==TRUE){
        ss<-seq(-10,-1,length=1000)
        ss<-10^ss
        n<-nrow(X)
        nn<-n- floor(n/k)
        lambda<-ss*nn*ncol(X)
    }
	p<-ncol(X)
	coefficients.jackknife=NULL
       if (compute.jackknife==TRUE){
   coefficients.jackknife<-array(dim=c(p,length(lambda),k))
	}
    
    n<-nrow(X)
if (is.null(groups) == FALSE) {
        f = as.factor(groups)
        k = length(levels(f))
        my.names = levels(f)
        all.folds <- split(1:n, f)
    }
    if (is.null(groups) == TRUE) {
        f <- rep(1:k, length = n)
        my.names <- 1:k
       all.folds <- split(sample(1:n), f)
    }
    cv.error.matrix<-cor.error.matrix<-matrix(0,k,length(lambda))
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain=X[-omit,,drop=FALSE]
        ytrain=y[-omit]
        Xtest=X[omit,,drop=FALSE]
        ytest=y[omit]
        ll<-lm.ridge(ytrain~Xtrain,scale=scale,lambda=lambda)
        coef.ll<-coef(ll)
	if (compute.jackknife==TRUE){
        	coefficients.jackknife[,,i]<-t(coef.ll[,-1])
	}        
	res<-matrix(,length(ytest),length(lambda))
        pred<-t(matrix(coef.ll[,1],nrow=length(lambda),ncol=length(ytest))) + Xtest%*%t(coef.ll[,-1])
	res<-pred-matrix(ytest,nrow=length(ytest),ncol=length(lambda))
	cv.error.matrix[i,]= apply(res^2,2,mean)  
        for (j in 1:length(lambda)){   
	cor.error.matrix[i,j]<-cor(pred[,j],ytest)
        }
    }
    cv.error<-apply(cv.error.matrix,2,mean)
    cor.error<-apply(cor.error.matrix,2,mean)
    lambda.opt<-lambda[which.min(cv.error)]
    lambda.opt.cor<-lambda[which.max(cor.error)]
    if (plot.it==TRUE){
        plot(lambda,cv.error,type="l",ylim="mean squared error")
    }
    rr<-lm.ridge(y~X,scale=scale,lambda=lambda.opt)
    coefficients<-coef(rr)
    intercept<-coefficients[1]
    coefficients<-coefficients[-1]
    rr.cor<-lm.ridge(y~X,scale=scale,lambda=lambda.opt.cor)
    coefficients.cor<-coef(rr.cor)
    intercept.cor<-coefficients.cor[1]
    coefficients.cor<-coefficients.cor[-1]
    return(list(cv.error=cv.error,cor.error=cor.error,cv.error.matrix=cv.error.matrix,cor.error.matrix=cor.error.matrix,intercept=intercept,coefficients=coefficients,lambda.opt=lambda.opt,intercept.cor=intercept.cor,coefficients.cor=coefficients.cor,lambda.opt.cor=lambda.opt.cor,coefficients.jackknife=coefficients.jackknife))
}
