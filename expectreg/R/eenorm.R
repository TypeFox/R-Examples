eenorm <-
function (y, main = "Normal E-E Plot", xlab = "Theoretical Expectiles", 
    ylab = "Sample Expectiles", plot.it = TRUE, datax = FALSE, 
    ...) 
{
    ey = expectile(y, seq(0.01, 0.99, by = 0.01))
    en = enorm(seq(0.01, 0.99, by = 0.01), mean(y), sd(y))
    if (plot.it) {
        if (datax) 
            plot(ey, en, main = main, xlab = ylab, ylab = xlab, 
                ...)
        else plot(en, ey, main = main, xlab = xlab, ylab = ylab, 
            ...)
        abline(0, 1, col = "red")
    }
    list(x = en, y = as.vector(ey))
}
