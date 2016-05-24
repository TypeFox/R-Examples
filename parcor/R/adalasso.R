adalasso<-function(X, y,k=10,use.Gram=TRUE,both=TRUE,intercept=TRUE){
    colnames(X)=1:ncol(X)
    n<-length(y)
    cv.adalasso<-NULL
    globalfit<-mylars(X,y,k=k,use.Gram=use.Gram,normalize=TRUE,intercept=intercept)
    coefficients.lasso=globalfit$coefficients
    intercept.lasso=globalfit$intercept
    cv.lasso<-globalfit$cv.lasso
    lambda<-globalfit$lambda
    lambda.lasso<-globalfit$lambda.opt
    coefficients.adalasso=NULL
    lambda.adalasso<-intercept.adalasso<-NULL
    if (use.Gram==TRUE){
        type="covariance"
    }
    if (use.Gram==FALSE){
        type="naive"
    }
    if (both==TRUE){ 
      # cross-validation for adaptive lasso
        all.folds <- split(sample(1:n),rep(1:k,length=n))
        residmat <- matrix(0, length(lambda), k)
    
        for (i in seq(k)) {
            omit <- all.folds[[i]]
            Xtrain<-X[-omit,,drop=FALSE]
            ytrain<-y[-omit]
            Xtest<-X[omit,,drop=FALSE]
            ytest<-y[omit]
            my.lars<-mylars(Xtrain,ytrain,k=k,normalize=TRUE,use.Gram=use.Gram,intercept=intercept)
            coef.lasso<-my.lars$coefficients
            weights <- 1/abs(coef.lasso[ abs(coef.lasso)>0 ])
            #cat(paste("-- non-zero weights ",length(weights),"\n"))
            if (length(weights)==0){
                residmat[,i]<-mean((mean(ytrain)-ytest)^2)
            }
            if (length(weights)==1){
                  residmat[,i]=mean((ytest -my.lars$intercept - Xtest%*%coef.lasso)^2)
            }
            if (length(weights)>1){
                XXtrain <- Xtrain[ , names(weights), drop=FALSE]
                XXtest<-Xtest[ , names(weights), drop=FALSE]
                XXtrain <- scale(XXtrain, center=FALSE, scale=weights)
                XXtest<-scale(XXtest, center=FALSE, scale=weights)
                #cat(paste("ncol of XXtrain: ",ncol(XXtrain),"\n"))
                fit<-glmnet(XXtrain,ytrain,type.gaussian=type,standardize=FALSE,intercept=intercept)
                pred<-predict(fit, newx=XXtest, type = "response",s = lambda)
                if (length(omit) == 1){
                    pred <- matrix(pred, nrow = 1)
                }
                residmat[, i] <- apply((ytest - pred)^2, 2, mean)
            }
    }
    cv <- apply(residmat, 1, mean)
    cv.adalasso<-min(cv)
    weights <- 1/abs(coefficients.lasso[ abs(coefficients.lasso)>0 ])
    coefficients.adalasso<-rep(0,ncol(X))
    names(coefficients.adalasso)<-1:ncol(X)
    if (length(weights)>0){
        XX <- X[ , names(weights), drop=FALSE]
        if ( length(weights)==1 )  XX <- XX/weights        
        else  XX <- scale(XX, center=FALSE, scale=weights)
        if (length(weights)<=1){
            intercept.adalasso=intercept.lasso 
            coefficients.adalasso<-coefficients.lasso
            lambda.adalasso=0
        }
        else{
    fit<-glmnet(XX,y,type.gaussian=type,standardize=FALSE,intercept=intercept)
    lambda.adalasso<-lambda[which.min(cv)]
    coefficients=predict(fit,type="coefficients",s=lambda.adalasso)
    intercept.adalasso<-coefficients[1]
    coefficients.adalasso[names(weights)]<-coefficients[-1]/weights
        }
    }
    }
    return(list(cv.lasso=cv.lasso,lambda.lasso=lambda.lasso,cv.adalasso=cv.adalasso,lambda.adalasso=lambda.adalasso,intercept.lasso=intercept.lasso, intercept.adalasso=intercept.adalasso, coefficients.lasso=coefficients.lasso,coefficients.adalasso=coefficients.adalasso))
}
