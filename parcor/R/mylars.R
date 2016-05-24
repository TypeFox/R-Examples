
mylars<-function (X, y, k = 10,use.Gram=TRUE,normalize=TRUE,intercept=TRUE) 
{
    x<-X
    n<-length(y)
    all.folds <- split(sample(1:n),rep(1:k,length=n))
    
    if (use.Gram==TRUE){
        type="covariance"
    }
    if (use.Gram==FALSE){
        type="naive"
    }
    globalfit<-glmnet(x,y,family="gaussian",standardize=normalize,type.gaussian=type,intercept=intercept)
    lambda<-globalfit$lambda
    residmat <- matrix(0, length(lambda), k)
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        fit <- glmnet(x[-omit, ,drop=FALSE], y[-omit],type.gaussian=type,standardize=normalize,family="gaussian",intercept=intercept)
        fit <- predict(fit, newx=x[omit, , drop = FALSE], type = "response", 
            s = lambda)
        if (length(omit) == 1) 
            fit <- matrix(fit, nrow = 1)
        residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
    }
    cv <- apply(residmat, 1, mean)
    cv.lasso<-min(cv)
    cv.error <- sqrt(apply(residmat, 1, var)/k)
    lambda.opt<-lambda[which.min(cv)]
    coefficients=predict(globalfit,type="coefficients",s=lambda.opt)
    inter=coefficients[1]
    coefficients=coefficients[-1]
    names(coefficients)=1:ncol(X)
    object <- list(lambda=lambda,cv=cv,lambda.opt=lambda.opt,cv.lasso=cv.lasso,intercept=inter,coefficients=coefficients)
    invisible(object)
}
