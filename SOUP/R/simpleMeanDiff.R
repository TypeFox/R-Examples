#-- meanDiff statistics in PP.R  --#
.simpleMeanDiff <- function(dataset, groups, indexMat)
{
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
    
    #- observed statistics
    Ttemp <- t(dataset) %*% CM
    T[1, , ] <- Ttemp

    #- permutation statistics
    for(bb in 2L:(B + 1))
    {
        ind <- indexMat[, bb - 1]
        data.p <- dataset[ind, , drop = FALSE]
        Ttemp <- crossprod(data.p, CM)
        T[bb, , ] <- Ttemp
    }# END:for-bb
    ##- last permutation==observed    
#    T[B + 1, , ] <- T[1, , ]

    dimnames(T) <- list(
        c("p-obs", paste("p-*", seq_len(B), sep = "")), colnames(dataset), labsPC
    )
    return(T)
}#=END=
