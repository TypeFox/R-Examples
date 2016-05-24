HWLratio <- function (X, verbose = TRUE, x.linked = FALSE) 
{
if(!x.linked) { # autosomal marker
    if (length(X) != 3 | any(X < 0)) 
        stop("HWLratio: X is not a 3 by 1 non-negative count vector")
    if (any(!is.wholenumber(X))) {
        warning("Genotype counts are not integers, counts will be rounded.")
        X <- round(X, digits = 0)
    }
    Xhom <- X[homozyg(X)]
    Xhet <- X[heterozyg(X)]
    X <- c(min(Xhom), Xhet, max(Xhom))
    nAA <- X[1]
    nAB <- X[2]
    nBB <- X[3]
    n <- sum(X)
    nA <- 2 * nAA + nAB
    nB <- 2 * nBB + nAB
    nxxvec <- X
    lnnxxvec <- log(X)
    ind <- nxxvec != 0
    nxxvec <- nxxvec[ind]
    lnnxxvec <- lnnxxvec[ind]
    LambdaA <- sum(nxxvec * lnnxxvec)
    nxvec <- c(nA, nB)
    lnnxvec <- log(nxvec)
    ind <- nxvec != 0
    nxvec <- nxvec[ind]
    lnnxvec <- lnnxvec[ind]
    LambdaB <- sum(nxvec * lnnxvec)
    lnLambda <- nAB * log(2) + LambdaB - n * log(4 * n) - LambdaA
    Lambda <- exp(lnLambda)
    G2 <- -2 * lnLambda
    pval <- pchisq(G2, 1, lower.tail = FALSE)
    if (verbose) {
        cat("Likelihood ratio test for Hardy-Weinberg equilibrium\n")
        cat("G2 =", G2, "DF = 1", "p-value =", pval, "\n")
    }
} else { # x-linked marker
    if (length(X) != 5 | any(X < 0)) 
        stop("HWLratio: X is not a 5 by 1 non-negative count vector for an x-linked marker")
    if (any(!is.wholenumber(X))) {
        warning("Genotype counts are not integers, counts will be rounded.")
        X <- round(X, digits = 0)
    }
    lab <- names(X)
    if(!all(lab %in% c("A","AA","AB","B","BB")))
        stop("Unknown genotypes occurred. Supply counts as a named vector c(A,AA,AB,B,BB)")
    n <- sum(X)         
    nfAA <- X[lab=="AA"]
    nfAB <- X[lab=="AB"]
    nfBB <- X[lab=="BB"]
    nmA <- X[lab=="A"]
    nmB <- X[lab=="B"]
    nm <- nmA+nmB
    nf <- n - nm
    X <- c(nmA,nmB,nfAA,nfAB,nfBB)
    nA <- nmA + 2*nfAA + nfAB
    nB <- nmB + 2*nfBB + nfAB
    nt <- nA+nB

    nxxvec <- X # genotype counts
    lnnxxvec <- log(X)
    ind <- nxxvec != 0
    nxxvec <- nxxvec[ind]
    lnnxxvec <- lnnxxvec[ind]
    LambdaA <- -sum(nxxvec * lnnxxvec)
         
    nxvec <- c(nA, nB) # allele counts
    lnnxvec <- log(nxvec)
    ind <- nxvec != 0
    nxvec <- nxvec[ind]
    lnnxvec <- lnnxvec[ind]
    LambdaB <- sum(nxvec * lnnxvec)
         
    nsexvec <- c(nm,nf) # sex counts
    lnnsexvec <- log(nsexvec)
    ind <- nsexvec != 0
    nsexvec <- nsexvec[ind]    
    lnnsexvec <- lnnsexvec[ind]    
    LambdaC <- sum(nsexvec * lnnsexvec)
        
    lnLambda <- LambdaA + LambdaB + LambdaC + nfAB*log(2) - nt*log(nt)
    DF <- 2
    
    Lambda <- exp(lnLambda)
    G2 <- -2 * lnLambda
    pval <- pchisq(G2, DF, lower.tail = FALSE)
    if (verbose) {
        cat("Likelihood ratio test for Hardy-Weinberg equilibrium for an X-linked marker\n")
        cat("G2 =", G2, "DF =", DF, "p-value =", pval, "\n")
    }         
}
    out <- list(Lambda = Lambda, G2 = G2, pval = pval)
}

