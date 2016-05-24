plot.csBiFunction <-
function (x, ...) 
{
    n <- length(x$alpha)
    plot(1:n, x$alpha, type = "l", ylab = "Functions", xlab = "Time", 
        ...)
    lines(1:n, x$beta, lty = 2)
}
