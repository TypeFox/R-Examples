#' Two-stage hybrid LASSO model.
#'
#' This function performs a LASSO logistic regression model using a two-stage hybrid procedure, namely the TSLasso logistic regression model, produces an optimal set of predictors and returns the robust estimations of coefficients of the selected predictors.
#'
#' @param x predictor matrix.
#' @param y response variable, a factor object with values of 0 and 1.
#' @param lambda.candidates the lambda candidates in the cv.lqa function, with the default values from 0.001 to 5 by=0.01.
#' @param kfold the number of folds of cross validation - default is 10. Although kfold can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is kfold=3.
#' @param seed seed for random sampling, with the default value 0123.
#' @export
#' @import glmnet
#' @import SiZer
#' @import lqa
#' @references
#' [1] Guo, P., Zeng, F., Hu, X., Zhang, D., Zhu, S., Deng, Y., Hao, Y. (2015). Improved Variable 
#' Selection Algorithm Using a LASSO-Type Penalty, with an Application to Assessing Hepatitis B 
#' Infection Relevant Factors in Community Residents. PLoS One, 27;10(7):e0134151.
#' Zou, H. (2006). The Adaptive Lasso And Its Oracle Properties. Journal of the American Statistical 
#' Association, 101(476), 1418:1429.
#' @examples
#' library(datasets)
#' head(iris)
#' X <- as.matrix(subset(iris,iris$Species!="virginica")[,-5])
#' Y <- as.numeric(ifelse(subset(iris,iris$Species!="virginica")[,5]=='versicolor', 0, 1))
#' # Fit a two-stage hybrid LASSO (TSLasso) logistic regression model.
#' # The parameters of lambda.candidates in the following example are set as small values to  
#' # reduce the running time, however the default values are proposed.
#' TSLasso.fit <- TSLasso(x=X, y=Y, lambda.candidates=list(seq(0.1, 1, by=0.05)), kfold=3, seed=0123)
#' # Variables selected by the TSLasso model.
#' TSLasso.fit$var.selected
#' # Coefficients of the selected variables.
#' TSLasso.fit$var.coef
TSLasso=function(x, y, lambda.candidates=list(seq(0.001, 5, by=0.01)), kfold=10, seed=0123){
    set.seed(seed)
    varx <- colnames(x)
    rowx <- nrow(x)
    n <- length(y)
    if (rowx!=n){
      stop("The number of rows in x is not equal to the length of y!")
                 }
    cvfit <- cv.glmnet(x=x, y=y, type.measure="deviance", nfolds=kfold, family="binomial")
    model.final <- cvfit$glmnet.fit
    nzero <- as.matrix(coef(model.final, s=cvfit$lambda.min))
    nzero <- as.matrix(nzero[rownames(nzero)!="(Intercept)",])
    var_Lasso <- names(nzero[nzero[,1]!=0,])
    cat("The LASSO algorithm finish!", "\n")
    if(!is.null(var_Lasso)){
    # Adaptive Lasso model
    x_Ada <- x[, var_Lasso]
    y_Ada <- y
    w <- (abs(as.vector(nzero[nzero[,1]!=0,])))^(-1)
    cvfit_Ada <- cv.lqa(y_Ada, x_Ada, lambda.candidates=lambda.candidates, family=binomial(), 
                        penalty.family=lasso, n.fold=3, loss.func="aic.loss", intercept=TRUE)
    lamd <- cvfit_Ada$lambda.opt
    fit_Ada <- lqa(y_Ada~x_Ada, family=binomial(), penalty=adaptive.lasso(lambda=lamd, al.weights=w))
    var_Ada1 <- round(fit_Ada$coef, 4)
    var_Ada2 <- names(round(fit_Ada$coef, 4))
    var_Ada3 <- var_Ada1[var_Ada2!="(Intercept)"]
    names(var_Ada3) <- var_Lasso
    var_Ada4 <- names(var_Ada3[var_Ada3!=0])
    cat("The adaptive LASSO algorithm finish!", "\n")
                            }
    else {
      print("No variables are selected by the LASSO algorithm!")
      var_Ada4 <-c()
          }
 var_Ada1 <- var_Ada1[-1]
 names(var_Ada1) <- var_Ada4
 var_Ada1 <- -var_Ada1
 Myresult <- list(var.selected=var_Ada4, var.coef=var_Ada1)
 Myresult
}