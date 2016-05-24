PlotFreq <-
function (path.mcmc, ipop, iloc, iall, printit = FALSE, path) 
{
    fileparam <- paste(path.mcmc, "parameters.txt", sep = "")
    param <- as.matrix(read.table(fileparam))
    nit <- as.numeric(param[param[, 1] == "nit", 3])
    thinning <- as.numeric(param[param[, 1] == "thinning", 3])
    filter.null.alleles <- as.logical(param[param[, 1] == "filter.null.alleles", 
        3])
    allele.numbers <- c(scan(paste(path.mcmc, "allele.numbers.geno2.txt", 
        sep = "")), scan(paste(path.mcmc, "allele.numbers.geno1.txt", 
        sep = "")), scan(paste(path.mcmc, "number.levels.ql.txt", 
        sep = "")))
    nalmax <- max(allele.numbers)
    filef <- paste(path.mcmc, "frequencies.txt", sep = "")
    f <- as.matrix(read.table(filef))
    npopmax <- ncol(f)
    nloc <- length(allele.numbers)
    sub1 <- rep(FALSE, times = iall - 1)
    sub2 <- TRUE
    sub3 <- rep(FALSE, times = nalmax - iall)
    sub <- c(sub1, sub2, sub3)
    sub1 <- rep(FALSE, (iloc - 1) * nalmax)
    sub2 <- sub
    sub3 <- rep(FALSE, (nloc - iloc) * nalmax)
    sub <- rep(c(sub1, sub2, sub3), times = nit/thinning)
    plot(f[sub, ipop], xlab = paste("Index of MCMC iteration", 
        " (x ", thinning, ")", sep = ""), ylab = paste("Frequency of allele", 
        iall, "at locus", iloc), type = "l", ylim = c(0, 1))
    if (iall < allele.numbers[iloc]) {
        title(main = paste("Allele frequencies of allele", iall, 
            "for locus", iloc, "in population", ipop))
    }
    if (filter.null.alleles & (iall == allele.numbers[iloc])) {
        title(main = paste("Frequencies of null alleles", "for locus", 
            iloc, "in population", ipop))
    }
    if (printit == TRUE) {
        postscript(file = paste(path, "freq.pop", ipop, ".loc", 
            iloc, ".ps", sep = ""))
        par(mfrow = c(ceiling(sqrt(allele.numbers[iloc])), ceiling(sqrt(allele.numbers[iloc]))))
        sub1 <- rep(FALSE, times = iall - 1)
        sub2 <- TRUE
        sub3 <- rep(FALSE, times = nalmax - iall)
        sub <- c(sub1, sub2, sub3)
        sub1 <- rep(FALSE, (iloc - 1) * nalmax)
        sub2 <- sub
        sub3 <- rep(FALSE, (nloc - iloc) * nalmax)
        sub <- rep(c(sub1, sub2, sub3), times = nit/thinning)
        plot(f[sub, ipop], xlab = paste("Index of MCMC iteration", 
            " (x ", thinning, ")", sep = ""), ylab = paste("Frequency of allele", 
            iall, "at locus", iloc), type = "l", ylim = c(0, 
            1))
        if (iall < allele.numbers[iloc]) {
            title(main = paste("Allele frequencies of allele", 
                iall, "for locus", iloc, "in population", ipop))
        }
        if (filter.null.alleles & (iall == allele.numbers[iloc])) {
            title(main = paste("Frequencies of null alleles", 
                "for locus", iloc, "in population", ipop))
        }
        dev.off()
    }
}
