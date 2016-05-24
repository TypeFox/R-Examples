#--- NPC statistic for RCB design PARAMETRIC CASE          ---#
#--- (nonparametric case is just like nonparametric npc.SP ---#
.strataTtest <- function(dataset, groups, strata, indexMat, linearInter = FALSE) {
    #- global variables
    groups <- as.factor(groups)
    strata <- as.factor(strata)
    B <- NCOL(indexMat)
    nObs <- NROW(dataset)
    p <- NCOL(dataset)
    C <- length(tab.groups <- table(groups))
    S <- length(tab.str <- table(strata))
    K <- C * (C - 1)/2

    ##- labels
    labsMat <- t(outer(levels(groups), levels(groups), FUN = paste, sep = "-"))
    labsPC <- labsMat[lower.tri(labsMat)]
    
    #- matrix of the p.values resulting from t.test
    T  <- array(NA, c(B + 1, p, K))
    ##- Contrasts Matrix for test statistics
    CM <- array(0, c(nObs, K))
    for(ss in seq_len(S)) {
        tab <- table(groups[strata == levels(strata)[ss]])
        CM[strata == levels(strata)[ss], ] <- .DesM(tab) / rep.int(tab, tab)
    }# END:for-ss

    ## - Contrast Matrix last steps: dividing by correct num.s
    ##  and sum up together columns of different strata
    sumMat <- kronecker(rep.int(1L, S), diag(K))
    
    ##- elements for the parametric statistics (C-sample t-test, equal var.)
    ##  pairwise sample size sum Matrix (1/n_i + 1/n_j) for the t-test calculation
    pairDiffMat <- abs(.DesM(rep.int(1L, C)))
    pairwiseSize <- matrix(t(pairDiffMat) %*% (1/tab.groups), nrow = p, ncol = K, byrow = TRUE)

    ##- Design Matrix for the sigma calculation
    X.lm <- model.matrix(
        model.frame("~ strata + groups", data = data.frame(groups, strata)),
        contrasts = list(groups = "contr.treatment", strata = "contr.treatment"),
        data = data.frame(groups, strata)
    )# END:X.lm

    ##- linear interaction with strata levels
    if(linearInter) {
        dataset <- dataset / matrix(.unfactor(strata), nrow = nObs, ncol = p)
    }# END:if-linearInter

    ##- observed statistics
    Ttemp <- t(dataset) %*% CM
    ##  sigma for the t-test statistic
    # sig <- apply(dataset, MARGIN = 2, FUN = .sigma, X.lm = X.lm) * sqrt(pairwiseSize)
    sig <- .sigma(dataset, X.lm = X.lm) * sqrt(pairwiseSize)
    Ptemp <- pt(Ttemp/sig, df = sum(tab.groups) - C - S, lower.tail = FALSE, log.p = FALSE)
    T[1, , ] <- Ptemp

    ##- permutation (bootstrap) statistics
    for(bb in 2L:(B + 1))
    {
        ind <- indexMat[, bb - 1, ]
        ##- permutations of rows of the dataset
        data.p <- dataset[ind[!is.na(ind)], , drop = FALSE]
        ##- test statistics
        Ttemp <- t(data.p) %*% CM
        # sig <- apply(data.p, 2, FUN = .sigma, X.lm = X.lm) * sqrt(pairwiseSize)
        sig <- .sigma(dataset, X.lm = X.lm) * sqrt(pairwiseSize)
        Ptemp <- pt(Ttemp/sig, df = sum(tab.groups) - C - S, lower.tail = FALSE, log.p = FALSE)
        #- NonParametric Combination
        T[bb, , ] <- Ptemp
    }# END:for-bb
    #- last permutation = observed
#    T[B + 1, , ] <- T[1, , ]
    dimnames(T) <- list(
        c("p-obs", paste("p-*", seq_len(B), sep = "")), colnames(dataset), labsPC
    )
    return(T)
}#=END=
