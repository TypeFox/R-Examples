mylasso<-function(x,y,lambda,standardize=TRUE,intercept=TRUE)
# compute the lasso estimator for a given tuning parameter lambda
{
    globalfit<-glmnet(x,y,standardize=standardize,intercept=intercept)
    fitlasso = predict(globalfit, type = "coefficients", s = lambda)
    beta0<-fitlasso[1]
    beta<-fitlasso[-1]
    object<-list()
    object$beta0<-beta0
    object$beta<-beta
    object$lambda<-lambda
    object
}