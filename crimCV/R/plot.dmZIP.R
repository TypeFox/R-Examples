plot.dmZIP <-
function (x, ...) 
{
    prob <- x$prob
    Xb <- x$X %*% x$beta
    Zg <- x$Z %*% x$gamma
    lambda <- exp(Xb)
    p <- exp(Zg)
    p <- p/(1 + p)
    mu <- (1 - p) * lambda
    tt <- 1:nrow(mu)
    matplot(tt, mu, type = "l", lty = 1, lwd = 2, xlab = "time", 
        ylab = "muhat(time)")
}
