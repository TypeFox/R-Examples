#---    NPC statistic for Secondary Performances    ---#
.strataMeanDiff <- function(dataset, groups, strata, indexMat, linearInter = FALSE) {
    
    ##- global variables
    groups <- as.factor(groups)
    strata <- as.factor(strata)
    B <- NCOL(indexMat)
    nObs <- NROW(dataset)
    p <- NCOL(dataset)
    C <- length(tab.groups <- table(groups))
    S <- length(tab.strata <- table(strata))
    K <- C * (C - 1)/2
    
    ##- labels
    labsMat <- t(outer(levels(groups), levels(groups), FUN = paste, sep = "-"))
    labsPC  <- labsMat[lower.tri(labsMat)]

    ##- matrices of test statistics
    T  <- array(NA, c(B + 1, p, K * S))
    #>
    # browser()
    #<
    ##- Contrasts Matrix for the test statistics, column positions inside CM
    #  are given by '((ss - 1)*K + 1):(ss * K)', row positions are given by
    #  'strata == levels(strata)[ss]'
    CM <- array(0, c(nObs, S * K))
    for(ss in seq_len(S)) {
        tab <- table(groups[strata == levels(strata)[ss]])
        CM[strata == levels(strata)[ss], ((ss - 1)*K + 1):(ss * K)] <-
            .DesM(tab) / rep.int(tab, tab)
    }# END:for-ss

    #- summation matrix for sum over strata
    sumMat <- kronecker(rep(1L, S), diag(K))
    
    ##- linear interaction with strata levels
    if(linearInter) {
        dataset <- dataset / matrix(.unfactor(strata), nrow = nObs, ncol = p)
    }# END:if-linearInter

    ##- observed statistics
    tmp <- t(dataset) %*% CM
    T[1, , ] <- tmp
    #>
    # browser()
    #<
    #- permutation statistics
    data.p <- dataset
    for(bb in 2L:(B + 1)) {
        ind <- indexMat[, bb - 1, ]
        #- permutations of rows of the dataset
        data.p <- dataset[ind[!is.na(ind)], , drop = FALSE]
        tmp <- crossprod(data.p, CM)
        T[bb, , ] <- tmp
    }# END:for-bb
    ##- last permutation = observed
#    T[B + 1, , ] <- T[1, , ]
    ##- sum over strata 
    T2 <- tensor(T, sumMat, 3, 1)
    dimnames(T2) <- list(
        c("p-obs", paste("p-*", seq_len(B), sep = "")), colnames(dataset), labsPC
    )
    return(T2)
}#=END=
