rgl <-
function (n, lambda1 = 0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
    param = "fmkl") 
{
    lambdas <- .gl.parameter.tidy(lambda1, lambda2, lambda3, 
        lambda4, param)
    if (!fun.check.gld(lambdas, param = param)) {
        return(NA)
    }
    p <- runif(n)
    res <- qgl(p, lambda1 = lambdas, param = param)
    res
}
