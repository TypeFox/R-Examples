##-----------------------------------------##
##  MODIFICARE IN MODO DA CALCOLARE SOLO   ##
##  I COEFFICIENTI E I RESIDUI DELL' 'lm'  ##
##-----------------------------------------##
##  NON FUNZIONA CON CONTRASTI DIVERSI DA  ##
##  QUELLI "TREATMENT" PER COSTRUZIONE     ##
##-----------------------------------------##
.strataLmCoef <- function(dataset, groups, strata, indexMat, deg, linearInter = FALSE,
    tY = function(x){x}, tX = function(x){x}) {
    #-----------------------------------------------------------------#
    ##- global variables
    groups <- as.factor(groups)
    B    <- NCOL(indexMat)# + 1
    nObs <- NROW(dataset)
    p    <- NCOL(dataset)
    C    <- length(tab.groups <- table(groups))
    S    <- length(tab.str <- table(strata))
    K    <- C * (C - 1)/2
    ##- Naming columns of dataset if absent
    if(is.null(colnames(dataset))) {
        colnames(dataset) <- LETTERS[1L:p]
    }# END:if-null-varNames
    ##- de-factorize strata, it must be numeric in the case of interaction
    strata <- .unfactor(strata)
    ##- checks missings
    if(missing(deg)) {
        deg <- S - 1
    }# END:if-missing-degreeOfPolynomial

    ##- labels
    labsMat <- t(outer(levels(groups), levels(groups), FUN = paste, sep = "-"))
    labsPC <- labsMat[lower.tri(labsMat)]

    ##- matrix of test statistics
    T <- array(NA, c(B + 1, p, K))
    ##- p.values matrix in the npc.PP style
    Ptemp <- array(NA, dim = c(p, K))

    ##- contrasts matrix for groups
    groupsContr <- model.matrix(
        model.frame("~ groups", data = data.frame(groups)),
        contrasts = list(groups = "contr.treatment"), data = data.frame(groups)
    )
    
    ##- Linear interaction between levels of strata and groups
    if(linearInter) {
        groupsContr[, -1] <- groupsContr[, -1] * matrix(strata, nrow = nObs, ncol = C - 1)
    }# END:if-inter
    
    ##- polynomials transformation of strata
    strataMat <- array(0, dim = c(nObs, deg))
    for(gg in seq_len(deg)) {
        strataMat[, gg] <- tX(strata)^gg
    }# END:for-gg
    if(deg > 1) {
        colnames(strataMat) <- c("strata", paste("strata.", 2L:deg, sep = ""))
    } else{ 
        colnames(strataMat) <- "strata"
    }# END:if-colnames    
    
    
    ##- matrix for pairwise differences of coefficients
    pairDiffCoef <- rbind(
        rep.int(0L, K),
        cbind(diag(C - 1), -.DesM(rep.int(1L, C - 1))),
        array(0, c(NCOL(strataMat), K))
    )
    
    ##- matrices for constrained beta    
    # Z <- cbind(groupsContr[, 1], strata.mat, groupsContr[, -1])
    Z <- cbind(groupsContr, strataMat)
    ZZ <- crossprod(Z)# %*% Z
    Z.inv <- solve(ZZ)
    Dtemp <- round(1/diag(crossprod(pairDiffCoef, Z.inv) %*% pairDiffCoef))
    ##- useful for calculates beta pairwise differences
    if(length(Dtemp) == 1) {
        Z.inv2 <- Z.inv %*% pairDiffCoef * Dtemp
    } else {
        Z.inv2 <- Z.inv %*% pairDiffCoef %*% diag(Dtemp)
    }# END:if-Dtemp
    
    ##- observed statistics
    ## applying transformation to response variables
    dataset <- apply(dataset, MARGIN = 2, FUN = tY)
    # X <- data.frame(apply(dataset, 2, FUN = tY), groupsContr, strataMat)    
    qrX <- qr(Z)
    F.df2 <- nrow(qrX$qr) - ncol(qrX$qr)

    ##  coefficients and residuals of the 'lm' estimation
    e2 <- diag(crossprod(qr.resid(qrX, dataset)))
    pairBeta <- crossprod(pairDiffCoef, qr.coef(qrX, dataset))
    signPairBeta <- t(pairBeta > 0)
    
    ##- Loop over variables
    for(i in seq_len(p)) {
        pBeta <- pairBeta[, i]
        ##  differences between constrained and unconstrained coefficients
        if(length(Dtemp) == 1) {
            betaDiff <- -(Z.inv2 * c(pBeta))
        } else {
            betaDiff <- -(Z.inv2 %*% diag(c(pBeta)))
        }# END:if-Dtemp
        ##  sum of squares in the constrained case
        em2 <- diag(crossprod(betaDiff, Z.inv) %*% betaDiff)
        ##- F test
        Ptemp[i, ] <- pt(
            sqrt(F.df2 * em2 / e2[i]),
            df = F.df2, lower.tail = FALSE, log.p = FALSE
        ) / 2
    }# END:for-i
    Ptemp[signPairBeta] <- 1 - Ptemp[signPairBeta]
    T[1, , ] <- Ptemp

    ##- permutation statistics
    for(bb in 2L:(B + 1))
    {
        ind <- indexMat[, bb - 1, ]
        data.p <- dataset[ind[!is.na(ind)], , drop = FALSE]
        ##  alternative 'lm' estimation (only coefs and residuals)
        e2 <- diag(crossprod(qr.resid(qrX, data.p)))
        pairBeta <- crossprod(pairDiffCoef, qr.coef(qrX, data.p))
        signPairBeta <- t(pairBeta > 0)
        ## Loop over variables
        for(i in seq_len(p)) {
            pBeta <- pairBeta[, i]
            if(length(Dtemp) == 1) {
                betaDiff <- -(Z.inv2 * c(pBeta))
            } else {
                betaDiff <- -(Z.inv2 %*% diag(c(pBeta)))
            }# END:if-Dtemp
            ##- differences "constr.beta - beta"
            # em2 <- e2[i] + diag(crossprod(betaDiff, ZZ) %*% betaDiff)
            em2 <- diag(crossprod(betaDiff, Z.inv) %*% betaDiff)
            ##- F test
            Ptemp[i, ] <- pt(
                sqrt(F.df2 * em2 / e2[i]),
                df = F.df2, lower.tail = FALSE, log.p = FALSE
            ) / 2
        }# END:for-i
        Ptemp[signPairBeta] <- 1 - Ptemp[signPairBeta]
        T[bb, , ] <- Ptemp
    }# END:for-bb
    ##- last permutation == observed
#    T[B + 1, , ] <- T[1, , ]
    dimnames(T) <- list(
        c("p-obs", paste("p-*", seq_len(B), sep = "")), colnames(dataset), labsPC
    )
    return(T)
}#=END=

