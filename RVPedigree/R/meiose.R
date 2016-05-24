meiose <- function(i, indices.fam, indices.haplotypes, haplotypes.perm)
{
    if (runif(1) > 0.5) {
        return(haplotypes.perm[indices.fam[i, 1], indices.haplotypes])
    } else {
        return(haplotypes.perm[indices.fam[i, 2], indices.haplotypes])
    }
}
