EstimateFreqNA <-
function (path.mcmc) 
{
    fileparam <- paste(path.mcmc, "parameters.txt", sep = "")
    param <- as.matrix(read.table(fileparam))
    nit <- as.numeric(param[param[, 1] == "nit", 3])
    thinning <- as.numeric(param[param[, 1] == "thinning", 3])
    filter.null.alleles <- as.logical(param[param[, 1] == "filter.null.alleles", 
        3])
    if (!filter.null.alleles) {
        stop("The null alleles option was not activated before MCMC inference")
    }
    allele.numbers.dip <- scan(paste(path.mcmc, "allele.numbers.geno2.txt", 
        sep = ""))
    allele.numbers.hap <- scan(paste(path.mcmc, "allele.numbers.geno1.txt", 
        sep = ""))
    nb.levels.ql <- scan(paste(path.mcmc, "number.levels.ql.txt", 
        sep = ""))
    filef <- paste(path.mcmc, "frequencies.txt", sep = "")
    f <- as.matrix(read.table(filef))
    npopmax <- ncol(f)
    nlocd <- length(allele.numbers.dip)
    nloch <- length(allele.numbers.hap)
    nql <- length(nb.levels.ql)
    ncolt <- nlocd + nloch + nql
    res <- matrix(nrow = nlocd, ncol = 2, data = 0)
    for (ipop in 1:npopmax) {
        for (iloc in 1:nlocd) {
            iall <- allele.numbers.dip[iloc]
            sub1 <- rep(FALSE, times = iall - 1)
            sub2 <- TRUE
            sub3 <- rep(FALSE, times = max(allele.numbers.dip) - 
                iall)
            sub <- c(sub1, sub2, sub3)
            sub1 <- rep(FALSE, (iloc - 1) * max(allele.numbers.dip))
            sub2 <- sub
            sub3 <- rep(FALSE, (ncolt - iloc) * max(allele.numbers.dip))
            sub <- rep(c(sub1, sub2, sub3), times = nit/thinning)
            ff <- f[sub, ipop]
            ff <- ff[ff > 0]
            res[iloc, 1] <- res[iloc, 1] + sum(ff)
            res[iloc, 2] <- res[iloc, 2] + length(ff)
        }
    }
    for (iloc in 1:nlocd) {
        res[iloc, 1] <- res[iloc, 1]/res[iloc, 2]
    }
    rownames(res) <- paste("Locus", 1:nlocd)
    res[, 1]
}
