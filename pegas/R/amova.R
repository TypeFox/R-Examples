## amova.R (2015-12-23)

##   Analysis of Molecular Variance

## Copyright 2010-2015 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

amova <- function(formula, data = NULL, nperm = 1000, is.squared = FALSE)
{
    y.nms <- as.character(as.expression(formula[[2]]))
    rhs <- formula[[3]]
    gr.nms <- as.character(as.expression(rhs))
    ## keep the highest level first:
    if (length(rhs) > 1) gr.nms <- unlist(strsplit(gr.nms, "/"))

    data.env <- if (is.null(data)) environment(formula)
                else as.environment(data)

    if (any(sapply(gr.nms, function(x) ! is.factor(get(x, envir = data.env)))))
        warning("elements in the rhs of the formula are not all factors")

    ## fix by Zhian Kamvar (2015-04-30)
    gr <- as.data.frame(sapply(gr.nms, get, envir = data.env))
    y <- get(y.nms, envir = environment(formula))
    if (any(is.na(y)))
        warning("at least one missing value in the distance object.")
    if (!is.squared) y <- y^2 # square the distances
    if (class(y) == "dist") y <- as.matrix(y)
    if (!is.matrix(y))
        stop("the lhs of the formula must be either a matrix or an object of class 'dist'.")
    n <- dim(y)[1] # number of individuals
    Nlv <- length(gr) # number of levels

### 5 local functions
    ## a simplified version of tapply(X, INDEX, FUN = sum):
    foo <- function(x, index) unlist(lapply(split(x, index), sum))

    getSSD <- function(y, gr, Nlv, N, n) {
        SSD <- numeric(Nlv + 2)
        SSD[Nlv + 2] <- sum(y/(2 * n)) # total SSD
        ## calculate SSD *within* each level:
        for (i in 1:Nlv) {
            p <- gr[, i] # extract the grouping at this level
            SSD[i + 1] <- sum((y/(2 * N[[i]])[p])[outer(p, p, "==")])
        }
        ## now differentiate to get the SSDs:
        if (Nlv > 1)
            for (i in 2:Nlv) SSD[i] <- SSD[i] - SSD[i + 1]
###          for (i in Nlv:2) SSD[i] <- SSD[i] - SSD[i + 1] # old line replaced by the previous one on 2015-12-22
        SSD[1] <- SSD[Nlv + 2] - sum(SSD[-(Nlv + 2)])
        SSD
    }

    getDF <- function(gr, Nlv, N, n) {
        df <- numeric(Nlv + 2)
        df[1:Nlv] <- unlist(lapply(N, length)) # number of groups
        ## get the df from the lowest level to the highest one:
        df[Nlv + 1] <- n
        for (i in (Nlv + 1):2) df[i] <-  df[i] - df[i - 1]
        df[1] <- df[1] - 1
        df[Nlv + 2] <- n - 1 # total df
        df
    }

    getNcoefficient <- function(gr, Nlv, N, n) {
        Nig <- N[[Nlv]] # numbers from the lowest level (ie, pop)
        Nig2 <- Nig^2
        npop <- length(Nig)
        if (Nlv == 1) # do classic variance components estimation
            ncoef <- (n - sum(Nig2)/n)/(npop - 1)
        else {
            if (Nlv == 2) { # use eq.9a-c in Excoffier et al. 1992
                ncoef <- numeric(3) # n, n', and n'', respectively
                G <- nlevels(gr[, 1]) # the number of groups
                ## use the fact that match() returns the 1st match of the
                ## 1st argument into the 2nd one, so this returns a factor
                ## indicating to which group each pop belongs to
                g <- gr[, 1][match(1:npop, as.integer(gr[, 2]))]
                npopBYgr <- tabulate(g) # number of pops in each group
                A <- sum(foo(Nig2, g)/foo(Nig, g))
                ncoef[1] <- (n - A)/sum(npopBYgr - 1) # corrects eq.9a in Excoffier et al. 1992
                ncoef[2] <- (A - sum(Nig2)/n)/(G - 1)
                ncoef[3] <- (n - sum(foo(Nig, g)^2/n))/(G - 1)
            } else { # use approximate formulae
                ncoef <- numeric(Nlv + 1)
                ## the coefficients are indexed from the highest to the lowest
                ## level to ease computation of the variance components
                ncoef[Nlv] <- (n - sum(Nig2)/n)/(npop - 1) # same than above
                ncoef[Nlv + 1] <- 1 # for the within pop level
                for (i in 1:(Nlv - 1)) {
                    group <- gr[, i]
                    g <- group[match(1:npop, as.integer(gr[, i + 1]))]
                    A <- sum(foo(Nig, g)^2)/sum(foo(Nig, g))
                    ncoef[i] <- (n - A)/(nlevels(group) - 1)
                }
            }
        }
        names(ncoef) <- letters[1:length(ncoef)] # suggestion by Rodney Dyer
        ncoef
    }

    getVarComp <- function(MSD, Nlv, ncoef) {
        ## get the variance components bottom-up
        if (Nlv == 1)
            sigma2 <- c((MSD[1] - MSD[2])/ncoef, MSD[2]) # 'n' changed to 'ncoef' (fix by Qixin He, 2010-07-15)
        else {
            sigma2 <- numeric(Nlv + 1)
            if (Nlv == 2) {
                sigma2[3] <- MSD[3]
                sigma2[2] <- (MSD[2] - sigma2[3])/ncoef[1]
                sigma2[1] <- (MSD[1] - MSD[3] - ncoef[2]*sigma2[2])/ncoef[3]
            } else {
                sigma2[Nlv + 1] <- MSD[Nlv + 1]
                for (i in Nlv:1) {
                    sel <- i:(Nlv + 1)
                    sigma2[i] <- (MSD[i] - sum(ncoef[sel]*sigma2[sel]))/ncoef[i]
                }
            }
        }
        names(sigma2) <- c(names(gr), "Error")
        sigma2
    }

    ## fit the AMOVA model to the data:
    N <- lapply(gr, tabulate) # number of individuals in each group
    SSD <- getSSD(y, gr, Nlv, N, n)
    df <- getDF(gr, Nlv, N, n)
    MSD <- SSD/df
    ncoef <- getNcoefficient(gr, Nlv, N, n)
    sigma2 <- getVarComp(MSD, Nlv, ncoef)

    ## output the results:
    res <- list(tab = data.frame(SSD = SSD, MSD = MSD, df = df,
                row.names = c(names(gr), "Error", "Total")),
                varcoef = ncoef, varcomp = sigma2, call = match.call())
    class(res) <- "amova"

### Below, "pop" is used for the lowest level of the geo structure,
### and "region" the one just above pop.
    if (nperm) {
        rSigma2 <- matrix(0, nperm, length(sigma2)) # Note that length(sigma2) == Nlv + 1
        ## First shuffle all individuals among all pops to assess the
        ## within-pop var components; if there is only one level, this
        ## is used to test both var components and no need to go further.
        j <- if (Nlv == 1) 1:2 else Nlv + 1
        for (i in 1:nperm) {
            rY <- perm.rowscols(y, n)
            ## the hierarchical structure is not changed so
            ## need to recalculate only the SSD and var comp
            rSSD <- getSSD(rY, gr, Nlv, N, n)
            rSigma2[i, j] <- getVarComp(rSSD/df, Nlv, ncoef)[j]
        }
        if (Nlv > 1) {
            ## for the lowest level (i.e., pop), we just permute individuals within the
            ## level just above, so the hierarchical structure is unchanged
            j <- Nlv
            ## 'L' contains the indices of the individuals for each region
            L <- lapply(levels(gr[, j - 1]), function(x) which(gr[, j - 1] == x))
            for (i in 1:nperm) {
                rind <- unlist(lapply(L, sample))
                rY <- y[rind, rind]
                rSSD <- getSSD(rY, gr, Nlv, N, n)
                rSigma2[i, j] <- getVarComp(rSSD/df, Nlv, ncoef)[j]
            }
            if (Nlv > 2) {
                ## code used from the 2nd lowest level up to the penultimate one
                for (j in (Nlv - 1):2) {
                    above <- gr[, j - 1]
                    L <- lapply(levels(above), function(x) which(above == x))
                    for (i in 1:nperm) {
                        rind <- integer(0)
                        for (k in L) rind <- c(rind, sample(k))
                        rind <- unlist(lapply(L, sample))
                        rY <- y[rind, rind]
                        rGR <- gr[rind, ]
                        rN <- lapply(rGR, tabulate)
                        rSSD <- getSSD(rY, rGR, Nlv, rN, n)
                        rDF <- getDF(rGR, Nlv, rN, n)
                        rNcoef <- getNcoefficient(rGR, Nlv, rN, n)
                        rSigma2[i, j] <- getVarComp(rSSD/rDF, Nlv, rNcoef)[j]
                    }
                }
            }
            ## for the highest level permute the level below
            N2 <- N[[2]]
            Higr <- gr[, 1][cumsum(N2)] # find at what highest level the 2nd-most one belongs to
            rGR <- gr # copy gr whose 1st column will be changed below
            for (i in 1:nperm) {
                ## with mapply(rep, sample....) the structure below the highest level is kept
                ## and the structures are assigned from the 2nd level within the 1st one
                rGR[, 1] <- unlist(mapply(rep, sample(Higr), each = N2, SIMPLIFY = FALSE))
                rN <- lapply(rGR, tabulate)
                rSSD <- getSSD(rY, rGR, Nlv, rN, n)
                rDF <- getDF(rGR, Nlv, rN, n)
                rNcoef <- getNcoefficient(rGR, Nlv, rN, n)
                rSigma2[i, 1] <- getVarComp(rSSD/rDF, Nlv, rNcoef)[1]
            }
        }
        P <- numeric(Nlv + 1)
        for (j in 1:(Nlv + 1))
            P[j] <- sum(rSigma2[, j] >= sigma2[j])/(nperm + 1)
        P[Nlv + 1] <- NA
        res$varcomp <- data.frame(sigma2 = res$varcomp, P.value = P)
    }
    res
}

print.amova <- function(x, ...)
{
    cat("\n\tAnalysis of Molecular Variance\n\nCall: ")
    print(x$call)
    cat("\n")
    print(x$tab)
    cat("\nVariance components:\n")
    if (is.data.frame(x$varcomp)) {
        x$varcomp["Error", "P.value"] <- NA
        printCoefmat(x$varcomp, na.print = "")
    } else print(x$varcomp)
    cat("\nVariance coefficients:\n")
    print(x$varcoef)
    cat("\n")
}
