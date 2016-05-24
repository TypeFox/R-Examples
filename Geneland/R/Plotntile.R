Plotntile <-
function (path.mcmc, burnin, printit = FALSE, file) 
{
    fileparam <- paste(path.mcmc, "parameters.txt", sep = "/")
    param <- as.matrix(read.table(fileparam))
    thinning <- as.numeric(param[param[, 1] == "thinning", 3])
    filenpp <- paste(path.mcmc, "nuclei.numbers.txt", sep = "")
    npp <- scan(filenpp)
    if (burnin > 0) {
        sub <- -(1:burnin)
    }
    else {
        sub <- 1:length(npp)
    }
    if (printit == TRUE) {
        postscript(file)
        par(mfrow = c(1, 2))
        plot((1:length(npp)) * thinning, npp, type = "l", ylab = "Number of Tiles", 
            xlab = paste("Index of MCMC iteration", " \n Whole chain", 
                sep = ""))
        hist(npp[sub], plot = TRUE, prob = TRUE, xlab = paste("Number of tiles along the chain \n(after a burnin of ", 
            burnin, "x", thinning, " it.)", sep = ""), main = "Number of tiles along the chain  \n  after burnin")
        dev.off()
    }
    else {
        par(mfrow = c(1, 2))
        plot((1:length(npp)) * thinning, npp, type = "l", ylab = "Number of Tiles", 
            xlab = paste("Index of MCMC iteration", "\n Whole chain", 
                sep = ""))
        hist(npp[sub], plot = TRUE, prob = TRUE, xlab = paste("Number of tiles along the chain \n(after a burnin of ", 
            burnin, "x", thinning, " it.)", sep = ""), main = "Number of tiles along the chain  \n  after burnin")
    }
}
