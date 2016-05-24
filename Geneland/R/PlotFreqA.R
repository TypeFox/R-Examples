PlotFreqA <-
function (path.mcmc, iloc, iall, printit = FALSE, path) 
{
    fileparam <- paste(path.mcmc, "parameters.txt", sep = "/")
    param <- as.matrix(read.table(fileparam))
    thinning <- as.numeric(param[param[, 1] == "thinning", 3])
    filter.null.alleles <- as.logical(param[param[, 1] == "filter.null.alleles", 
        3])
    allele.numbers <- c(scan(paste(path.mcmc, "allele.numbers.geno2.txt", 
        sep = "")), scan(paste(path.mcmc, "allele.numbers.geno1.txt", 
        sep = "")), scan(paste(path.mcmc, "number.levels.ql.txt", 
        sep = "")))
    nalmax <- max(allele.numbers)
    nloc <- length(allele.numbers)
    filefa <- paste(path.mcmc, "ancestral.frequencies.txt", sep = "")
    fa <- as.matrix(read.table(filefa))
    dev.new()
    plot(fa[, (iloc - 1) * nalmax + iall], xlab = paste("Index of MCMC iteration", 
        " (x ", thinning, ")", sep = ""), ylab = paste("Frequency of allele", 
        iall, "at locus", iloc), type = "l", ylim = c(0, 1))
    title(main = ifelse(iall == 1, paste("Allele freq. in ancestral pop. at locus", 
        iloc), ""))
    if (printit == TRUE) {
        postscript(file = paste(path, "freq.ancestral.pop.loc", 
            iloc, ".ps", sep = ""))
        plot(fa[, (iloc - 1) * nalmax + iall], xlab = paste("Index of MCMC iteration", 
            " (x ", thinning, ")", sep = ""), ylab = paste("Frequency of allele", 
            iall, "at locus", iloc), type = "l", ylim = c(0, 
            1))
        title(main = ifelse(iall == 1, paste("Allele freq. in ancestral pop. at locus", 
            iloc), ""))
        dev.off()
    }
}
