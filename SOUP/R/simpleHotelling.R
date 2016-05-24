.simpleHotelling <- function(dataset, groups, indexMat, p.valuesType) {
    #-- Hotelling statistics in simple.aov --#
    #-- BALANCED EXPERIMENTS ONLY
    ##- Local utility functions
    ##  quadratic form
    .quadForm <- function(x, A) {
        return(crossprod(x, A) %*% x)
    }# END:.quadForm_equivalent to diag(crossprod(t(crossprod(muDiff, S)), muDiff))
    .cholScale <- function(x, A) {
        return(
            sum(tcrossprod(x, chol(A)))
        )
    }# END:.cholScale_equivalent to rowSums(tcrossprod(t(muDiff), chol(A)))
    ##- Global variables
    dataset.orig <- dataset
    groups <- as.factor(groups)
    B <- NCOL(indexMat)
    p <- NCOL(dataset)
    nObs <- NROW(dataset)
    C <- length(tab <- table(groups))
    n <- min(tab)
    K <- C * (C - 1)/2

    ##- labels
    labsMat <- t(outer(levels(groups), levels(groups), FUN = paste, sep = "-"))
    labsPC <- labsMat[lower.tri(labsMat)]

    ##- checking unbalancing of the experiment
    balanced <- all(tab == n)
    if(!balanced) {
        groups.ub <- groups
        nObs.ub <- nObs
        ind.ub <- unlist(tapply(1:nObs.ub, groups.ub,
            FUN = .mySample, size = n))
        dataset <- dataset.orig[ind.ub,]
        nObs <- NROW(dataset)
        groups <- groups.ub[ind.ub]
        tab <- table(groups)
    }# END:if!balanced

    ##- covariance matrix
    S <- matrix(0, nrow = p, ncol = p)
    Sg  <- array(NA, dim = c(p, p, C))
    
    ##- signs matrix
    sgn <- array(NA, dim = c(B + 1, K))
    ##- matrix of raw test-statisics
    T <- array(NA, dim = c(B + 1, K))
    
    #- Contrasts Matrix for the pairwise differences of means
    CM <- .DesM(tab) / rep(tab, tab)
    
    ##- correction factor for attaining the asymptotic F distribution is:
    ##  balanced case   = ((n^2)/(2 * n * p)) * (n*C - p - 1)/(n*C - C)
    ##  UNbalanced case = ((n_i * n_j)/(p*(n_i + n_j))) * (nObs - p - 1)/(nObs - C)
    corrFactor <- (n/(2 * p)) * (n*C - p - 1)/(n*C - C)
    #
    #=START=
    #- matrix of means
    muDiff <- corrFactor * (t(dataset) %*% CM)
    #- variance-covariance matrix
    # for(cc in seq_len(C)) {
        # Sg[, , cc] <- (n - 1) * cov(
            # x = dataset[groups == levels(groups)[cc], ], y = NULL, na.method = 4L, FALSE
        # )
    # }# END:for
    # S <- solve(apply(Sg, MARGIN = c(1, 2), FUN = sum) / (n*C - C))
    for(cc in seq_len(C)) {
        S <- S + (n - 1) * cov(x = dataset[groups == levels(groups)[cc], ], 
            y = NULL)
    }# END:for
    S <- solve(S/(n*C - C))
    #<--- INSERIRE QUI CALCOLO n_i * n_j PER CASO SBILANCIATO
    #- observed statistics
    # sgn <- sign(apply(muDiff, 2, FUN = .cholScale, A = S))
    sgn[1, ] <- sign(rowSums(tcrossprod(t(muDiff), chol(S))))
    # Stats <- apply(corrFactor * muDiff, 2, FUN = .quadForm, A = S)
    T[1, ] <- diag(crossprod(crossprod(S, muDiff), muDiff))
    #>
    # browser()
    #<
    #- permutation (bootstrap) statistics
    data.p <- dataset
    for(bb in 2L:(B + 1))
    {
        ##- rebalancing
        if(!balanced) {
            ind.ub <- unlist(tapply(1L:nObs.ub, groups.ub, FUN = .mySample, size = n))
            dataset <- dataset.orig[ind.ub, , drop = FALSE]
        }# END:if
        ind <- indexMat[, bb - 1]
        data.p <- dataset[ind, , drop = FALSE]
        muDiff <- crossprod(data.p, CM)
        # for(cc in seq_len(C)) {
            # Sg[, , cc] <- (n - 1) * cov(
                # x = dataset[groups == levels(groups)[cc], ], y = NULL, na.method = 4L, FALSE
            # )
        # }# END:for
        # S <- solve(apply(Sg, MARGIN = c(1, 2), FUN = sum) / (n*C - C))
        ##- refresh and re-calculate covariance matrix
        S[] <- 0
        for(cc in seq_len(C))
        {
            S <- S + (n - 1) * cov(x = dataset[groups == levels(groups)[cc], ], 
                y = NULL)
        }# END:for
        S <- solve(S/(n*C - C))
        ##- sign of the test statistics using differences vector rescaled with cholesky
        ##  decomposition of S (inverse of covariance matrix)
        # sgn <- sign(apply(corrFactor * muDiff, 2, FUN = .cholScale, A = S))
        sgn[bb, ] <- sign(rowSums(tcrossprod(t(muDiff), chol(S))))
        # Stats <- apply(corrFactor * muDiff, 2, FUN = .quadForm, A = S)
        T[bb, ] <- diag(crossprod(t(crossprod(muDiff, S)), muDiff))
    }# END:for-bb
    #- last permutation = observed
#    sgn[B + 1, ] <- sgn[1, ]
#    T[B + 1, ] <- T[1, ]
    
    ##- mapping to asymptotic p-values or not
    if(p.valuesType == "asymptotic") {
        T <- pf(T, df1 = p, df2 = n*C - p - 1, lower.tail = FALSE, log.p = FALSE)/2
        T[sgn > 0] <- 1 - T[sgn > 0]
        ##- coherence check
        # pd <- matrix(NA, C, C)
        # pd[lower.tri(pd)] <- tmp
        # pd <- t(pd)
        # pd[lower.tri(pd)] <- 1 - tmp
        # ifelse(pd < .05, 1, 0)
    } else {
        T <- T * sgn
    }# END:ifelse-p.valuesType == "asymptotic"

    dimnames(T) <- list(
        c("p-obs", paste("p-*", seq_len(B), sep = "")), labsPC
    )
    return(T)
}#=END=
