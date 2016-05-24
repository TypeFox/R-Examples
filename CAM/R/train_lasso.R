train_lasso <-
function(X,y,pars = list())
{
    cvres <- cv.glmnet(X,y)
    mod <- glmnet(X,y,lambda = cvres$lambda.1se)
    result <- list()
    result$Yfit <- predict(mod,X)
    result$residuals <- y - result$Yfit
    result$model <- mod
    return(result)
}
