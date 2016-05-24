HWPerm <- function (x, nperm = 17000, verbose = TRUE, x.linked = FALSE, FUN = ifelse(x.linked,Chisquare.x,Chisquare), ...) 
{
if(!x.linked) { # autosomal marker
    n <- sum(x)
    nA <- 2 * x[1] + x[2]
    nB <- 2 * n - nA
    stat.obs <- FUN(x)
    pseudodist <- numeric(nperm)
    i1 <- seq(1, 2 * n, 2)
    i2 <- seq(2, 2 * n, 2)
    for (i in 1:nperm) {
        xx <- sample(c(rep("A", nA), rep("B", nB)))
        A1 <- xx[i1]
        A2 <- xx[i2]
        Geno <- paste(A1, A2, sep = "")
        Geno[Geno == "BA"] <- "AB"
        nAA <- sum(Geno == "AA")
        nAB <- sum(Geno == "AB")
        nBB <- sum(Geno == "BB")
        y <- c(AA = nAA, AB = nAB, BB = nBB)
        stat.pseudo <- FUN(y)
        pseudodist[i] <- stat.pseudo
    }
    nlarger <- sum(pseudodist >= stat.obs)
    pval <- nlarger/nperm
    if (verbose) {
        cat("Permutation test for Hardy-Weinberg equilibrium\n")
        cat("Observed statistic:", stat.obs, " ", nperm, "permutations. p-value:", 
            pval, "\n")
    }
} else { # X-linked marker
  if (length(x) != 5 | any(x < 0)) 
        stop("HWPerm: x is not a 5 by 1 non-negative count vector for an x-linked marker")
    if (any(!is.wholenumber(x))) {
        warning("Genotype counts are not integers, counts will be rounded.")
        x <- round(x, digits = 0)
    }
    lab <- names(x)
    if(!all(lab %in% c("A","AA","AB","B","BB")))
        stop("Unknown genotypes occurred. Supply counts as a named vector like c(A,AA,AB,B,BB)")
    n <- sum(x)         
    nfAA <- x[lab=="AA"]
    nfAB <- x[lab=="AB"]
    nfBB <- x[lab=="BB"]
    nmA <- x[lab=="A"]
    nmB <- x[lab=="B"]
    nm <- nmA+nmB
    nf <- n - nm
    x <- c(nmA,nmB,nfAA,nfAB,nfBB)
    nA <- nmA + 2*nfAA + nfAB
    nB <- nmB + 2*nfBB + nfAB
    nt <- nA+nB
    stat.obs <- FUN(x)
    pseudodist <- numeric(nperm)
    for (i in 1:nperm) {
        xx <- sample(c(rep("A", nA), rep("B", nB)))
        males <- xx[1:nm]
        nmAsim <- sum(males=="A")
        nmBsim <- sum(males=="B")
        females <- xx[(nm+1):nt]
        i1 <- seq(1, 2 * nf, 2)
        i2 <- seq(2, 2 * nf, 2)        
        A1 <- females[i1]
        A2 <- females[i2]
        Geno <- paste(A1, A2, sep = "")
        Geno[Geno == "BA"] <- "AB"
        nfAAsim <- sum(Geno == "AA")
        nfABsim <- sum(Geno == "AB")
        nfBBsim <- sum(Geno == "BB")
        y <- c(A = nmAsim, B = nmBsim, AA = nfAAsim, AB = nfABsim, BB = nfBBsim)
        stat.pseudo <- FUN(y)
        pseudodist[i] <- stat.pseudo
    }
    nlarger <- sum(pseudodist >= stat.obs)
    pval <- nlarger/nperm
    if (verbose) {
        cat("Permutation test for Hardy-Weinberg equilibrium of an X-linked marker\n")
        cat("Observed statistic:", stat.obs, " ", nperm, "permutations. p-value:", 
            pval, "\n")
    }         
}
    out <- list(stat = stat.obs, pval = pval)
}
