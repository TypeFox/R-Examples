"jackknife"<-
function(x, theta, ...)
{
    call <- match.call()
    n <- length(x)
    u <- rep(0, n)
    for(i in 1:n) {
        u[i] <- theta(x[ - i], ...)
    }
    thetahat <- theta(x, ...)
    jack.bias <- (n - 1) * (mean(u) - thetahat)
    jack.se <- sqrt(((n - 1)/n) * sum((u - mean(u))^2))
    return(list(jack.se=jack.se, 
                jack.bias=jack.bias, 
                jack.values = u, 
                call=call))
}
