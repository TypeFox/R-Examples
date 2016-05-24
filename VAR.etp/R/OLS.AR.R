OLS.AR <-
function(x,p)
{
    x <- as.matrix(x)
    n <- nrow(x)
    B <- LSM(x,p)
    b <- B$coef;xmat=B$xmat
    e <- B$resid    
    s2 <- sum(e^2)/(nrow(x)-length(b))

return(list(coef=b,resid=e,covmat=B$covmat,xmat=xmat))
}
