
mice.impute.norm3 <- function (y, ry, x, ridge = 10^(-5) , ...) 
{
    x <- cbind(1, as.matrix(x))
    parm <- .norm.draw3(y, ry, x, ridge=ridge ,  ...)
    return(x[!ry, ] %*% parm$beta + stats::rnorm(sum(!ry)) * parm$sigma)
}