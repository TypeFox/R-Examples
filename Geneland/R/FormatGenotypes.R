FormatGenotypes <-
function (genotypes, ploidy) 
{
    formatted <- genotypes
    nall <- numeric(ncol(genotypes)/2)
    if (ploidy == 2) {
        for (ic in 1:(ncol(genotypes)/2)) {
            ens <- sort(unique(c(genotypes[, 2 * ic - 1], genotypes[, 
                2 * ic])))
            max.ens <- length(ens)
            for (il in 1:(dim(genotypes)[1])) {
                formatted[il, 2 * ic - 1] <- ifelse(is.na(genotypes[il, 
                  2 * ic - 1]), NA, (1:max.ens)[genotypes[il, 
                  2 * ic - 1] == ens])
                formatted[il, 2 * ic] <- ifelse(is.na(genotypes[il, 
                  2 * ic]), NA, (1:max.ens)[genotypes[il, 2 * 
                  ic] == ens])
            }
            nall[ic] <- max.ens
        }
    }
    if (ploidy == 1) {
        for (ic in 1:(ncol(genotypes))) {
            ens <- sort(unique(genotypes[, ic]))
            max.ens <- length(ens)
            for (il in 1:(dim(genotypes)[1])) {
                formatted[il, ic] <- ifelse(is.na(genotypes[il, 
                  ic]), NA, (1:max.ens)[genotypes[il, ic] == 
                  ens])
            }
            nall[ic] <- max.ens
        }
    }
    formatted <- as.matrix(formatted)
    list(genotypes = formatted, allele.numbers = nall)
}
