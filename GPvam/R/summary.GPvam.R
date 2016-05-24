summary.GPvam <-
function (object, ...) 
{
    alpha <- 0.05
    c.level <- qnorm(1 - alpha/2)
    nyear <- object$num.year
    k <- dim(object$parameters)[1]
    AIC <- 2 * k - 2 * object$loglik
    AICc <- AIC + 2 * k * (k + 1)/(object$num.obs - k - 1)
    res2 <- c(list( AIC = AIC, AICc = AICc), as.list(object))
    class(res2) <- c("summary.GPvam")
    return(res2)
}
