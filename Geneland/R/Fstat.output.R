Fstat.output <-
function (coordinates = NULL, genotypes, ploidy = 2, burnin = NULL, 
    path.mcmc) 
{
    param <- as.matrix(read.table(paste(path.mcmc, "parameters.txt", 
        sep = "")))
    post.param <- as.matrix(read.table(paste(path.mcmc, "postprocess.parameters.txt", 
        sep = "")))
    pmpi <- read.table(paste(path.mcmc, "modal.pop.indiv.txt", 
        sep = ""))[, 3]
    npop <- scan(paste(path.mcmc, "populations.numbers.txt", 
        sep = ""))
    burnin <- as.numeric(post.param[post.param[, 1] == "burnin", 
        3])
    npop.max <- as.numeric(param[param[, 1] == "npopmax", 3])
    npop.est <- order(hist(npop[-(1:burnin)], breaks = seq(0.5, 
        npop.max + 0.5, 1), plot = FALSE)$counts, decreasing = TRUE)[1]
    Fstat(genotypes, npop.est, pmpi, ploidy = ploidy)
}
