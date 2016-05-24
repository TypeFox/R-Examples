PlotDrift <-
function (path.mcmc, printit = FALSE, file) 
{
    fileparam <- paste(path.mcmc, "parameters.txt", sep = "/")
    param <- as.matrix(read.table(fileparam))
    thinning <- as.numeric(param[param[, 1] == "thinning", 3])
    filedrift <- paste(path.mcmc, "drifts.txt", sep = "")
    drift <- as.matrix(read.table(filedrift))
    npopmax <- ncol(drift)
    if (printit == TRUE) {
        postscript(file = file)
        par(mfrow = c(ceiling(sqrt(npopmax)), ceiling(sqrt(npopmax))))
        for (iclass in 1:npopmax) {
            plot(drift[, iclass], xlab = paste("Index of MCMC iteration", 
                " (x ", thinning, ")", sep = ""), ylab = paste("Drift of population", 
                iclass), type = "l")
        }
        dev.off()
    }
    dev.new()
    par(mfrow = c(ceiling(sqrt(npopmax)), ceiling(sqrt(npopmax))))
    for (iclass in 1:npopmax) {
        plot(drift[, iclass], xlab = paste("Index of MCMC iteration", 
            " (x ", thinning, ")", sep = ""), ylab = paste("Drift of population", 
            iclass), type = "l", ylim = c(0, max(drift)))
    }
}
