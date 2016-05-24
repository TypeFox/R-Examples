HWExactPower <- function (X, alternative = "two.sided", pvaluetype = "dost", 
    verbose = FALSE, theta = 4) 
{
  # computes the power of an exact test for HWE for the given data vector, p-value
  # end effect size (theta).
    n <- sum(X)
    Xhom <- X[homozyg(X)]
    Xhet <- X[heterozyg(X)]
    nA <- 2 * Xhom[1] + Xhet # allele counts
    nB <- 2 * n - nA
    MaxHet <- min(nA, nB) # maximum number of heterozygotes.
    if (MaxHet < 2) {
        pval <- 1
        prob <- 1
        pofthesample <- 1
        ind <- 1
    }
    else {
        ind <- match(X[2], seq(MaxHet%%2, MaxHet, 2))
#        enAB <- nA * nB/(2 * n - 1)
#        enAB <- round(enAB, digits = 0)
        pA <- nA/(2*n)
        EstimatedF <- ThetatoF(pA,theta)
        enAB <- as.integer(round(2*pA*(1-pA)*(1-EstimatedF)*n))
        if ((enAB%%2) != (MaxHet%%2)) 
            enAB <- enAB + 1
        nAA <- (nA - enAB)/2
        nBB <- (nB - enAB)/2
        initialprob <- 1
        AboveExp <- NULL
        BelowExp <- NULL
        if (enAB < MaxHet) 
            AboveExp <- CompProbUpPower(nAA, nBB, enAB, initialprob, MaxHet, theta)
        BelowExp <- CompProbDownPower(nAA, nBB, enAB, initialprob, theta)
        prob <- c(rev(BelowExp), initialprob, AboveExp)
        prob <- prob/sum(prob)
    }
    Plow <- cumsum(prob)
    Phigh <- 1 - c(0, Plow)
    Phigh <- Phigh[-length(Phigh)]
    Phwe <- pmin(1, 2 * Phigh, 2 * Plow)
    pofthesample <- prob[ind]
    pval <- switch(alternative,
                   greater = Phigh[ind],
                   less = Plow[ind], 
                   two.sided = switch(pvaluetype,
                                      dost = Phwe[ind],
                                      selome = sum(prob[prob <= pofthesample]),
                                      midp = sum(prob[prob < pofthesample])+0.5*pofthesample ,
                                      stop("invalid value for parameter pvaluetype")), 
                   stop("invalid value for parameter alternative"))
    if (verbose) {
        D <- HWChisq(X)$D
        cat("Haldane's Exact test for Hardy-Weinberg equilibrium\n")
        stringpvalue <- switch(pvaluetype, dost = "using DOST p-value\n", 
            selome = "using SELOME p-value\n", midp = "using MID p-value\n", stop("invalid value for parameter pvaluetype"))
        cat(stringpvalue)
        cat(paste("sample counts: n",names(Xhom[1])," = ",sep=""), Xhom[1], paste("n",names(Xhet)," = ",sep=""), Xhet,
            paste("n",names(Xhom[2])," = ",sep=""), Xhom[2], "\n")
        stringtwosided <- paste("H0: HWE (D==0), H1: D <> 0 \nD = ", 
            format(D, scientific = FALSE), "p = ", format(pval, 
                scientific = FALSE), "\n")  
        stringgreater <- paste("H0: HWE (D==0), H1: D > 0 \nD = ", 
            format(D, scientific = FALSE), "p = ", format(pval, 
                scientific = FALSE), "\n")
        stringless <- paste("H0: HWE (D==0), H1: D < 0 \nD = ", 
            format(D, scientific = FALSE), "p = ", format(pval, 
                scientific = FALSE), "\n")
        toprint <- switch(alternative, two.sided = stringtwosided, 
            greater = stringgreater, less = stringless)
        cat(toprint)
    }
    return(list(pval = pval, prob = prob, pofthesample = pofthesample))
}
