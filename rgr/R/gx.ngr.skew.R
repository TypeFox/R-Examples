gx.ngr.skew <-
function (xx) 
{
    x <- na.omit(xx)
    n <- length(x)
    xbar <- mean(x)
    xdif <- x - xbar
    skew <- sum(xdif^3) / sum(xdif^2)^1.5
    skew <- skew * ((n - 1)/n)^1.5
    invisible(skew)
}
