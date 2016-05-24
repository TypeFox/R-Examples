#-- Anderson-Darling statistics in SP.R  --#
.strataAD <- function(dataset, groups, strata, indexMat) {
    ##- global variables
    groups <- as.factor(groups)
    strata <- as.factor(strata)
    B <- NCOL(indexMat)
    nObs <- NROW(dataset)
    p <- NCOL(dataset)
    C <- length(table(groups))
    S <- length(table(strata))
    K <- C * (C - 1)/2

    ##- labels
    labsMat <- t(outer(levels(groups), levels(groups), FUN = paste, sep = "-"))
    labsPC  <- labsMat[lower.tri(labsMat)]

    ##- matrixes of the p.values resulting from t.test, results of NPC and
    ##   dataset re-allocated for pairwise comparisons
    T  <- array(NA, c(B + 1, p, K * S))
    T2 <- array(NA, c(B + 1, p, K))

    ##- Dummy variables Matrix, i.e. new dataset
    DM <- NULL
    for(j in seq_len(p)) {
        DM <- cbind(DM, .dummyze(dataset[,j]))
    }# END:for-dummyzing

    ##- number of categories of each variable
    Q <- as.integer(apply(dataset, MARGIN = 2, FUN = function(x) nlevels(factor(x))))

    #- Checking (un)balance of the experiment, Contrasts Matrix (CM) and Summation Matrix (SM)
    tab <- vector("list", S)
    CM <- array(0, c(nObs, S * K))
    SM <- array(0, c(nObs, S * K))
    for(ss in seq_len(S)) {
        tab <- table(groups[strata == levels(strata)[ss]])
        CM[strata == levels(strata)[ss], ((ss - 1)*K + 1):(ss * K)] <- 
            .DesM(tab) / rep.int(tab, tab)
        SM[strata == levels(strata)[ss], ((ss - 1)*K + 1):(ss * K)] <- 1/sum(tab)
    }# END:for-ss
    ##- change of sign: X > Y (in distributions) => F_x(t) < F_y(t) forall 't'
    CM <- -CM
    
    #- matrix for sums for sum over strata
    sumMat <- kronecker(rep(1L, S), diag(K))
    
    ##- Cumulative Sums matrix, block diagonal (CS.bd), and final Summation Matrix (fSM)
    CS <- diag(max(Q))
    CS[upper.tri(CS)] <- 1
    sq1 <- c(0L, cumsum(Q))
    CS.bd <- array(0, c(sum(Q), sum(Q)))
    fSM <- array(0, c(sum(Q) - length(Q), length(Q)))
    dd <- diag(p)
    for(k in seq_along(Q)) {
        CS.bd[(sq1[k] + 1):sq1[k + 1], (sq1[k] + 1):sq1[k + 1]] <- CS[1:Q[k], 1:Q[k]]
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
        ##- permutations of rows of DM
        ind <- indexMat[, bb - 1, ]
        FM.p <- FM[ind[!is.na(ind)], ]
        ## numerator (the last category is removed)
        num <- crossprod(CM, FM.p)[, -cumsum(Q)]
        ## denominator (the last category is removed)
        temp <- crossprod(SM, FM.p)[, -cumsum(Q)]
        den <- temp * (1 - temp)
        den[den < .Machine$double.eps] <- 0
        ## statistic
        stat <- t(num/sqrt(den))
        stat[!is.finite(stat)] <- 0
        T[bb, , ] <- t(crossprod(stat, fSM))
    }# END:for-bb
    ##- last permutation equal to the observed statistics
#    T[B + 1, , ] <- T[1, , ]
    
    #- sum over strata 
    T2[, , ] <- tensor(T, sumMat, 3, 1)
    dimnames(T2) <- list(
        c("p-obs", paste("p-*", seq_len(B), sep = "")), colnames(dataset), labsPC
    )
    return(T2)
}#=END=
