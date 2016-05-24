## Fst.R (2010-06-09)

##   F-Statistics

## Copyright 2009-2010 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

Fst <- function(x, pop = NULL)
{
    if (any(getPloidy(x) != 2))
        stop("Fst() requires diploid data")

    if (is.null(pop)) {
        pop <- x$population
        if (is.null(pop)) stop("no 'population' column in x")
    } else {
        pop <- if (is.numeric(pop)) x[, pop] else factor(pop)
    }
    r <- length(attr(pop, "levels")) # nb of pops
    pop <- as.integer(pop)
    nloci <- length(attr(x, "locicol"))
    nBYpop <- tabulate(pop)
    N <- nrow(x)
    nbar <- N/r

    ## get the frequencies in each population for each locus
    ALLELES <- getAlleles(x)
    p <- vector("list", nloci)
    for (j in 1:nloci)
        p[[j]] <- matrix(0, r, length(ALLELES[[j]]))
    h <- p
    for (i in 1:r) {
        s <- summary(x[pop == i, ]) # levels are preserved
        for (j in 1:nloci) {
            tmp <- s[[j]]
            p[[j]][i, ] <- tmp$allele
            allel <- names(tmp$allele)
            genot <- names(tmp$genotype)
            for (k in seq_along(allel)) {
                for (l in seq_along(genot)) {
                    ag <- unlist(strsplit(genot[l], "/"))
                    if (sum(ag %in% allel[k]) == 1)
                        h[[j]][i, k] <- h[[j]][i, k] + tmp$genotype[l]
                }
            }
        }
    }
    ## 'p' is a list with, for each locus, a matrix
    ##    with alleles as columns and populations as rows,
    ##    and its entries are the counts
    ## 'h' is the same but with the number of heterozygotes

    nC <- (N - sum(nBYpop^2)/N)/(r - 1)

    obj <- matrix(0, nloci, 3)
    for (j in 1:nloci) {
        ptild <- p[[j]]/(2 * nBYpop)
        pbar <- colSums(p[[j]])/(2 * N) # for each allele in the locus
        s2 <- colSums(nBYpop * (ptild - rep(pbar, each = r))^2)/((r - 1) * nbar)
        hbar <- colSums(h[[j]])/N # id.
        A <- pbar * (1 - pbar) - (r - 1) * s2/r
        a <- nbar * (s2 - (A - hbar/4)/(nbar - 1))/nC
        b <- nbar * (A - (2*nbar - 1) * hbar/(4*nbar))/(nbar - 1)
        c <- hbar/2
        obj[j, 1] <- 1 - sum(c)/sum(a + b + c)
        obj[j, 2] <- sum(a)/sum(a + b + c)
        obj[j, 3] <- 1 - sum(c)/sum(b + c)
    }
    dimnames(obj) <- list(names(x)[attr(x, "locicol")],
                          c("Fit", "Fst", "Fis"))
    obj
}
