## nuc.div.R (2015-10-27)

##   Nucleotide Diversity

## Copyright 2009-2015 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

nuc.div <- function(x, ...) UseMethod("nuc.div")

nuc.div.DNAbin <- function(x, variance = FALSE, pairwise.deletion = FALSE, ...)
{
    if (pairwise.deletion && variance) {
        warning("cannot compute the variance of nucleotidic diversity\nwith pairwise deletion: try 'pairwise.deletion = FALSE' instead.")
        variance <- FALSE
    }
    if (is.list(x)) x <- as.matrix(x)
    d <- dim(x)
    n <- d[1]
    ans <- sum(dist.dna(x, "RAW", pairwise.deletion = pairwise.deletion))/
        (n*(n - 1)/2)
    if (variance) {
        var <- (n + 1)*ans/(3*(n - 1)*d[2]) + 2*(n^2 + n + 3)*ans^2/(9*n*(n - 1))
        ans <- c(ans, var)
    }
    ans
}

nuc.div.haplotype <- function(x, variance = FALSE, pairwise.deletion = FALSE, ...)
{
    if (pairwise.deletion && variance) {
        warning("cannot compute the variance of nucleotidic diversity\nwith pairwise deletion: try 'pairwise.deletion = FALSE' instead.")
        variance <- FALSE
    }
    D <- dist.dna(x, "RAW", pairwise.deletion = pairwise.deletion)
    f <- sapply(attr(x, "index"), length)
    n <- sum(f)
    ff <- outer(f, f)
    ff <- ff[lower.tri(ff)]
    ans <- 2 * sum(ff * D) / (n * (n - 1)) # scale with n(n - 1) to have the same results than with nuc.div.DNAbin()
    if (variance) {
        var <- (n + 1)*ans/(3*(n - 1)*dim(x)[2]) + 2*(n^2 + n + 3)*ans^2/(9*n*(n - 1))
        ans <- c(ans, var)
    }
    ans
}
