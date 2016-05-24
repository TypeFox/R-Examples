train_linear <-
function(X,y,pars = list())
{
    mod <- lm(y ~ X)
    result <- list()
    result$Yfit = as.matrix(mod$fitted.values)
    result$residuals = as.matrix(mod$residuals)
    result$model = mod
    #for coefficients see list(mod$coef)
    return(result)
}
