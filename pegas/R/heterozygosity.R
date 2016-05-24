## heterozygosity.R (2002-08-28)

##   Heterozygosity at a Locus Using Gene Frequencies

## Copyright 2002 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

H <- function(x, variance = FALSE)
{
    if (!is.factor(x)) {
        if (is.numeric(x)) {
            n <- sum(x)
            k <- length(x)
            freq <- x/n
        } else x <- factor(x)
    }
    if (is.factor(x)) { # ne pas remplacer par `else'...
        n <- length(x)
        k <- nlevels(x)
        freq <- table(x)/n
    }
    sp2 <- sum(freq^2)
    H <- n * (1 - sp2) / (n - 1)
    if (variance) {
        sp3 <- sum(freq^3)
        var.H <- 2 * (2 * (n - 2) * (sp3 - sp2^2) + sp2 - sp2^2) / (n * (n - 1))
        return(c(H, var.H))
    }
    else return(H)
}

heterozygosity <- function(x, variance = FALSE) H(x, variance)
