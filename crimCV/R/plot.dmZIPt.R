plot.dmZIPt <-
function (x, ...) 
{
    prob <- x$prob
    Xb <- x$X %*% x$beta
    lambda <- exp(Xb)
    p <- exp(-x$tau * t(Xb))
    p <- t(p)
    p <- p/(1 + p)
    mu <- (1 - p) * lambda
    tt <- 1:nrow(mu)
    matplot(tt, mu, type = "l", lty = 1, lwd = 2, xlab = "time", 
        ylab = "muhat(time)")
}
