#-- t-test statistics in PP.R  --#
.simpleTtest <- function(dataset, groups, indexMat) {
    #- global variables
    groups <- as.factor(groups)
    B <- NCOL(indexMat)
    nObs <- NROW(dataset)
    p <- NCOL(dataset)
    C <- length(tab <- table(groups))
    K <- C * (C - 1)/2
    
    ##- labels
    labsMat <- t(outer(levels(groups), levels(groups), FUN = paste, sep = "-"))
    labsPC <- labsMat[lower.tri(labsMat)]
    
    #- matrix of statistics
    T <- array(NA, c(B + 1, p, K))
    
    #- Contrasts Matrix (CM), and checking (un)balance of the experiment
    CM <- .DesM(tab)/rep(tab, tab)
    
    ##- elements for the parametric statistics (C-sample t-test, equal var.)
    pairDiffMat <- abs(.DesM(rep.int(1L, C)))
    #- pairwise sample size sum Matrix (1/n_i + 1/n_j) for the t-test calculation
    pairwiseSize <- matrix(t(pairDiffMat) %*% (1/tab), nrow = p, ncol = K, byrow = TRUE)
    X.lm <- model.matrix(
        model.frame("~ groups", data = data.frame(groups)),
        contrasts = list(groups = "contr.treatment"), data = data.frame(groups)
    )
    #>
    # browser()
    #<
    #- observed statistics
    Ttemp <- t(dataset) %*% CM
    # sig <- apply(dataset, 2, FUN = .sigma, X.lm = X.lm) * sqrt(pairwiseSize)
    sig <- .sigma(dataset, X.lm = X.lm) * sqrt(pairwiseSize)
    Ptemp <- pt(Ttemp/sig, df = sum(tab) - C, lower.tail = FALSE, log.p = FALSE)
    T[1, , ] <- Ptemp

    #- permutation statistics
    for(bb in 2L:(B + 1))
    {
        ind <- indexMat[, bb - 1]
        data.p <- dataset[ind, , drop = FALSE]
        Ttemp <- t(data.p) %*% CM
        # sig <- apply(data.p, 2, FUN = .sigma, X.lm = X.lm) * sqrt(pairwiseSize)
        sig <- .sigma(dataset, X.lm = X.lm) * sqrt(pairwiseSize)
        Ptemp <- pt(Ttemp/sig, df = sum(tab) - C, lower.tail = FALSE, log.p = FALSE)
        T[bb, , ] <- Ptemp
    }# END:for-bb
    ##- last permutation==observed    
#    T[B + 1, , ] <- T[1, , ]

    dimnames(T) <- list(
        c("p-obs", paste("p-*", seq_len(B), sep = "")), colnames(dataset), labsPC
    )
    return(T)
}#=END=
