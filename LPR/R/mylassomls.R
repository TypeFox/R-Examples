mylassomls<-function(x,y,tau=0,lambda,standardize=TRUE,intercept=TRUE)
# compute Lasso+mLS estimator; lambda is the tuning parameter of the Lasso
{
    x<-as.matrix(x)
    n<-nrow(x)
    p<-ncol(x)
    globalfit<-glmnet(x,y,standardize=standardize,intercept=intercept)
    fitlasso <- predict(globalfit, type = "coefficients", s = lambda)
    betalasso<-fitlasso[-1]
    selectvar<-betalasso!=0
    beta0<-0
    beta<-rep(0,p)
    object<-list()
    if(sum(selectvar)>0){
    	ls.obj<-mls(x[,selectvar,drop=FALSE],y,tau,standardize,intercept)
    	beta0<-ls.obj$beta0
    	beta[selectvar]<-ls.obj$beta
    }
    object$beta0<-beta0
    object$beta<-beta
    object$lambda<-lambda
    object
}