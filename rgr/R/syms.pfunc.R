syms.pfunc <-
function () 
{
    frame()
    p <- c(0.2, 0.4, 0.6, 0.8, 1, 1.25, 1.6667, 2.5, 5)
    pl <- c("0.2", "0.4", "0.6", "0.8", "1", "1.25", "1.67", 
        "2.5", "5")
    np <- length(p)
    x <- 1:1000
    x <- x/1000
    y <- x
    plot(x, y, xlim = c(0, 1), xlab = "x", ylim = c(0, 1), ylab = "x^p", 
        type = "n", main = "Behaviour of Symbol Plotting Power Function")
    for (i in 1:np) {
        y <- x^p[i]
        lines(x, y, type = "l")
        text(x[500], y[500], pl[i], col = 2, cex = 1.5)
    }
    text(0, 1, "Function curves annotated by values of p", col = 2, 
        adj = 0)
    invisible()
}
