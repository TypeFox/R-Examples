## hw.test.R (2015-03-30)

##   Test of Hardy--Weinberg Equilibrium

## Copyright 2009-2014 Emmanuel Paradis, 2015 Thibaut Jombart

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

## This function outputs the multinomial coefficients with input:
## x: output from expand.genotype(, matrix = TRUE), see code in simulate.genotype()
.multinomialCoef <- function(x)
{
    n <- nrow(x)
    numerator <- factorial(ncol(x))
    coef <- numeric(n)
    for (i in seq_len(n))
        coef[i] <- numerator/prod(factorial(tabulate(x[i, ])))
    coef
}

proba.genotype <- function(alleles = c("1", "2"), p, ploidy = 2)
{
    o <- .sort.alleles(alleles, index.only = TRUE)
    alleles <- alleles[o]
    p <- if (missing(p)) rep(1/length(alleles), length(alleles)) else p[o]
    geno <- expand.genotype(length(alleles), ploidy = ploidy, matrix = TRUE)
    P <- .multinomialCoef(geno) * apply(geno, 1, function(i) prod(p[i]))
    names(P) <- apply(matrix(alleles[geno], ncol = ploidy), 1, paste, collapse = "/")
    P
}

expand.genotype <- function(n, alleles = NULL, ploidy = 2, matrix = FALSE)
{
    if (!is.null(alleles)) {
        alleles <- .sort.alleles(alleles)
        n <- length(alleles)
    }
    ans <- matrix(0L, 0L, ploidy)
    foo <- function(i, a) {
        for (x in a:n) {
            g[i] <<- x
            if (i < ploidy) foo(i + 1L, x) else ans <<- rbind(ans, g)
        }
    }
    g <- integer(ploidy)
    foo(1L, 1L)
    dimnames(ans) <- NULL
    if (is.character(alleles))
        ans <- matrix(alleles[ans], ncol = ploidy)
    if (!matrix)
        ans <- apply(ans, 1, paste, collapse = "/")
    ans
}


## CODE MODIFICATION BY THIBAUT JOMBART
## t.jombart@imperial.ac.uk
## March 2015
##

## generic
hw.test <- function (x, B=1000, ...){
    UseMethod("hw.test")
}

## method for class 'loci'
hw.test.loci <- function(x, B = 1000, ...)
{
    test.polyploid <- function(x, ploidy) {
        if (ploidy < 2) return(rep(NA_real_, 3))
        nms <- names(x$allele)
        all.prop <- x$allele/sum(x$allele)
        E <- proba.genotype(nms, all.prop, ploidy = ploidy)
        O <- E
        O[] <- 0
        O[names(x$genotype)] <- x$genotype
        E <- sum(O) * E
        chi2 <- sum((O - E)^2/E)
        DF <- length(E) - length(nms)
        c(chi2, DF, 1 - pchisq(chi2, DF))
    }
    y <- summary.loci(x)
    ploidy <- getPloidy(x)
    ans <- t(mapply(test.polyploid, y, ploidy = ploidy))
    dimnames(ans) <- list(names(y), c("chi^2", "df", "Pr(chi^2 >)"))
    if (B) {
        LN2 <- log(2)
        test.mc <- function(x, ploidy) {
            ## accept only diploids for the moment
            if (ploidy != 2) {
                warning("no Monte Carlo test available for polyploids")
                return(NA_real_)
            }
            n <- sum(x$genotype) # Nb of individuals
            Nall <- length(x$allele) # Nb of alleles
            all.geno <- expand.genotype(Nall)
            Ngeno <- length(all.geno)
            p <- numeric(B)
            ## the constant term below is dropped:
            ## lfactorial(n) - lfactorial(2*n) + sum(lfactorial(x$allele))
            ma <- unlist(mapply(rep, 1:Nall, x$allele))
            homoz.nms <- paste(1:Nall, 1:Nall, sep = "/") # to help find the homozygotes (see below)
            for (k in 1:B) {
                m <- sample(ma)
                dim(m) <- c(n, 2)
                ## check the order of both alleles:
                o <- m[, 1] > m[, 2]
                m[o, 1:2] <- m[o, 2:1]
                s <- factor(paste(m[, 1], m[, 2], sep = "/"), levels = all.geno)
                f <- tabulate(s, Ngeno)
                ## find the homozygotes:
                p[k] <- -sum(lfactorial(f)) + LN2 * sum(f[-match(homoz.nms, all.geno)])
            }
            ## find the homozygotes:
            h <- unlist(lapply(strsplit(names(x$genotype), "/"),
                               function(x) length(unique(x)))) == 1
            p0 <- -sum(lfactorial(x$genotype)) + LN2*sum(x$genotype[!h])
            sum(p <= p0)/B
        }
        ans <- cbind(ans, "Pr.exact" = mapply(test.mc, y, ploidy))
    }
    ans
} # end hw.test.loci



## method for class 'genind'
hw.test.genind <- function(x, B=1000, ...){
    ## checks
    byPop <- FALSE
    pop <- NULL
    if(is.null(pop)) pop <- pop(x)
    if(byPop && is.null(pop(x))){
        warning("byPop requested but pop(x) is NULL")
        byPop <- FALSE
    }

    ## convert genind to loci
    x <- as.loci(x)

    ## call hw.test.loci
    out <- hw.test.loci(x=x, B=B, ...)

    ## return result
    return(out)
} # end hw.test.genind




## version avec rhyper(): (Huber et al. 2006, Biometrics)
## (doesn't work)
#for (k in 1:B) {
#    N <- n
#    f <- x$allele
#    ff <- matrix(0, m, m)
#    for (i in 1:m) {
#        a <- rhyper(1, N, N, f[i]) #HG(2 * N, N, f[i])
#        y <- rhyper(1, a, N - a, f[i] - a)
#        ff[i, i] <- y #HG(N, a, f[i] - a)
#        N <- N - (f[i] - y)
#        f[i] <- f[i] - 2*y
#        b <- 2*N - f[i]
#        if (i != m) {
#            for (j in (i + 1):m) {
#                ff[i, j] <- rhyper(1, f[i], b - f[i], f[j]) #HG(b, f[i], f[j])
#                b <- b - f[j]
#                f[j] <- f[j] - ff[i, j]
#                f[i] <- f[i] - ff[i, j]
#            }
#        }
#    }
#    p[k] <- tmp - sum(lfactorial(ff)) + sum(diag(ff)) * LN2
#}
