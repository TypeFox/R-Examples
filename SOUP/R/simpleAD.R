.simpleAD <- function(dataset, groups, indexMat) {
    #-- Anderson-Darling statistics in PP.R  --#

    ##- global variables
    groups <- as.factor(groups)
    B <- NCOL(indexMat)
    nObs <- NROW(dataset)
    p <- NCOL(dataset)
    C <- length(tab <- table(groups))
    K <- C * (C - 1)/2

    ##- labels
    labsMat <- t(outer(levels(groups), levels(groups), FUN = paste, sep = "-"))
    labsPC <- labsMat[lower.tri(labsMat)]

    ##- matrix of statistics
    T <- array(NA, c(B + 1, p, K))

    ##- Dummy variables Matrix, i.e. new dataset
    DM <- NULL
    for(j in seq_len(p)) {
        DM <- cbind(DM, .dummyze(dataset[,j]))
    }# END:for-dummyze

    ##- number of categories of each variable
    Q <- as.integer(apply(dataset, MARGIN = 2, FUN = function(x) nlevels(factor(x))))
    
    ##- Contrasts Matrix (CM) and Summation Matrix (SM), and
    ##  checking (un)balance of the experiment
    CM <- -.DesM(tab) / rep(tab, tab)
    SM <- array(1/nObs, c(nObs, K))

    ##- Cumulative Sums matrix, block diagonal (CS.bd), and final Summation Matrix (fSM)
    CS <- diag(max(Q))
    CS[upper.tri(CS)] <- 1
    sq1 <- c(0L, cumsum(Q))
    CS.bd <- array(0, c(sum(Q), sum(Q)))
    fSM <- array(0, c(sum(Q) - length(Q), length(Q)))
    # sq2 <- c(0L, cumsum(Q - 1))
    dd <- diag(p)
    for(k in seq_along(Q)) {
        CS.bd[(sq1[k] + 1):sq1[k + 1], (sq1[k] + 1):sq1[k + 1]] <- CS[1:Q[k], 1:Q[k]]
        # fSM[(sq2[k] + 1):sq2[k + 1], k] <- 1
        fSM[, k] <- rep(dd[k, ], times = (Q - 1))
    }# END:for-Q

    ##- matrix of distribution function ('F') of each observation
    FM <- DM %*% CS.bd
    colnames(FM) <- colnames(DM)

    ##- observed statistics
    ## numerator (the last category is removed)
    num <- crossprod(CM, FM)[, -cumsum(Q)]
    ## denominator (the last category is removed)
    temp <- crossprod(SM, FM)[, -cumsum(Q)]
    den <- temp * (1 - temp)
    den[den < .Machine$double.eps] <- 0
    ## statistic
    stat <- num / sqrt(den)
    stat[!is.finite(stat)] <- 0
    T[1, , ] <- t(stat %*% fSM)
    
    ##- permutation statistics (it slightly differ from the
    ##  observed one for computational speed reasons)
    for(bb in 2L:(B + 1))
    {
        ind <- indexMat[, bb - 1]
        FM.p <- FM[ind, ]
        ## numerator
        num <- crossprod(CM, FM.p)[, -cumsum(Q)]
        ## denominator
        temp <- crossprod(SM, FM.p)[, -cumsum(Q)]
        den <- temp * (1 - temp)
        den[den < .Machine$double.eps] <- 0
        ## statistic
        stat <- t(num / sqrt(den))
        stat[!is.finite(stat)] <- 0
        T[bb, , ] <- t(crossprod(stat, fSM))
    }# END:for-bb
    ##- last permutation equal to the observed statistics
#    T[B + 1, , ] <- T[1, , ]
    dimnames(T) <- list(
        c("p-obs", paste("p-*", seq_len(B), sep = "")), colnames(dataset), labsPC
    )
    return(T)
}#=END=
