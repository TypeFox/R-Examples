plotBS <-
function (BS, alpha = 0.05, plot = TRUE, verbose = FALSE, main = "Bootstrap Histogram", 
    xlab = "Test Statistic Values", ylab = "Frequency") 
{
    B <- length(BS)
    if (plot == TRUE) {
        hist(BS[2:B], main = main, xlab = xlab, ylab = ylab)
        abline(v = BS[1])
    }
    if (verbose == TRUE) 
        cat("Realized Bootstrap is ", BS[1], "\n")
    p <- sum(BS[1] < BS[2:B])/B
    if (verbose == TRUE) {
        cat("p-value is ", p, "\n")
        if (p < alpha) 
            cat("Series was NOT stationary\n")
        else cat("Series was stationary\n")
    }
    return(p)
}
