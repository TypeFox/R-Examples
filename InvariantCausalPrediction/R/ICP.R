
ICP <- function (X, Y, ExpInd, alpha = 0.01, test = "normal", selection = c("lasso", 
    "all", "stability", "boosting")[if (ncol(X) <= 8) 2 else 4], 
    maxNoVariables = 8, maxNoVariablesSimult = 8, maxNoObs = 200, 
    showAcceptedSets = TRUE, showCompletion = TRUE, stopIfEmpty = FALSE, gof= max(0.01,alpha)) 
{
    if (is.vector(X) & is.numeric(X)) 
        X <- matrix(X, ncol = 1)
    if (!is.matrix(X) & !is.data.frame(X)) 
        stop("'X' must be a matrix or data frame")
    if (!is.vector(Y) & !is.factor(Y)) 
        stop("'Y' must be a vector or factor")
    if (is.function(test)) {
        pval <- test((1:10) + 0.5, 1:10)
        if (!is.numeric(pval)) 
            stop("function 'test' has to return a numeric value")
        if (length(pval) > 1) 
            stop("function 'test' needs to return a scalar (the p-value of the null hypothesis test that 'x' and 'z' are sampled from the same distribution")
        if (pval < 0 | pval > 1) 
            stop("the p-value of function 'test' needs to be in [0,1]")
    }
    if (!is.list(ExpInd)) {
        if (length(ExpInd) != length(Y)) 
            stop("if `ExpInd' is a vector, it needs to have the same length as `Y'")
        uni <- unique(ExpInd)
        if (length(uni) == 1) 
            stop(paste("there is just one environment ('ExpInd'=", 
                uni[1], " for all observations) and the method needs at least two distinct environments", 
                sep = ""))
        if (min(table(ExpInd)) <= 2) {
            cat("\n out put of 'table(ExpInd)':\n ")
            print(table(ExpInd))
            stop("one environment has just one or two observations (as supplied by 'ExpInd'); there need to be at least 3 (and ideally dozens) of observations in each environment; the output of 'table(ExpInd)' is given below to show the number of observations in each unique environment as supplied by 'ExpInd'")
        }
        K <- length(uni)
        ExpIndNEW <- list()
        for (uc in 1:K) {
            ExpIndNEW[[uc]] <- which(ExpInd == uni[uc])
            attr(ExpIndNEW[[uc]], "value") <- uni[uc]
        }
        ExpInd <- ExpIndNEW
        rm(ExpIndNEW)
    } else {
        ran <- range(unlist(ExpInd))
        if (ran[1] < 1) 
            stop(paste("if `ExpInd' is a list with indicies of observations, \n minimal entry has to be at least 1 but is", 
                ran[1]))
        if (ran[2] > length(Y)) 
            stop(paste("if `ExpInd' is a list with indicies of observations, \n maximal entry has to be at most equal \n to the length", 
                length(Y), "of the observations but is", ran[2]))
        if (min(sapply(ExpInd, length) <= 2)) 
            stop("one environment has just one or two observations (as supplied by 'ExpInd'); there need to be at least 3 (and ideally dozens) of observations in each environment")
    }
    getblanket <- getblanketall
    if (selection == "lasso") {
        getblanket <- getblanketlasso
    }
    if (selection == "stability") {
        getblanket <- getblanketstability
    }
    if (selection == "boosting") {
        getblanket <- getblanketboosting
    }
    if (is.data.frame(X)) {
        if (any(sapply(X, class) == "factor")) {
            Z <- X
            X <- matrix(nrow = nrow(Z), ncol = sum(sapply(X, 
                function(x) if (is.factor(x)) length(levels(x)) else 1)))
            cc <- 0
            colX <- character(0)
            for (k in 1:ncol(Z)) {
                if (is.numeric(Z[, k])) {
                  cc <- cc + 1
                  X[, cc] <- Z[, k]
                  colX <- c(colX, colnames(Z)[k])
                }
                else {
                  nf <- length(lev <- levels(Z[, k]))
                  for (nfc in 1:(nf - 1)) {
                    cc <- cc + 1
                    X[, cc] <- as.numeric(Z[, k] == lev[nfc])
                    colX <- c(colX, paste(colnames(Z)[k], "_", 
                      lev[nfc], sep = ""))
                  }
                }
            }
            X <- X[, 1:cc]
            colnames(X) <- colX
            X <- as.matrix(X)
        }
    }
    if (length(ucol <- unique(colnames(X))) < min(3, ncol(X))) 
        colnames(X) <- paste("Variable", 1:ncol(X), sep = "_")
    if (length(unique(Y)) == 2 & !is.factor(Y)) {
        warning("\n Y only has 2 unique values -- using classification")
        Y <- as.factor(Y)
    }
    K <- length(ExpInd)
    n <- nrow(X)
    p <- ncol(X)
    X <- cbind(rep(1, nrow(X)), X)
    ConfInt <- matrix(NA, nrow = 2, ncol = p)
    Coeff <- list()
    CoeffVar <- list()
    for (k in 1:p) {
        Coeff[[k]] <- numeric(0)
        CoeffVar[[k]] <- numeric(0)
    }
    pvall <- numeric(0)
    Pall <- numeric()
    cont <- TRUE
    pvalempty <- getpval(Y, X[, 1, drop = FALSE], ExpInd, test = test, 
        maxNoObs = maxNoObs)$pval
    if (pvalempty > alpha) {
        ConfInt <- matrix(0, nrow = 2, ncol = p)
        pvall <- c(pvall, pvalempty)
        if (showAcceptedSets) 
            cat(paste("\n accepted empty set"))
        if (stopIfEmpty) 
            cont <- FALSE
    }
    testsets <- if (any(!is.null(c(maxNoVariables, maxNoVariablesSimult)))){
        getblanket(X[, -1, drop = FALSE], Y, maxNoVariables = maxNoVariables, 
                   maxNoVariablesSimult = maxNoVariablesSimult)
    }else{
        getblanket(X[, -1, drop = FALSE], Y)
    }
    len <- sapply(testsets, length)
    lcsingle <- sum(len == 1)
    usedvariables <- unique(unlist(testsets))
    lc <- 0
    printoutat <- ifelse(length(testsets) > 0, 2^(1:ceiling(log2(length(testsets)))), 
        NA)
    if(cont & stopIfEmpty) intersection <- NULL
    while (cont && lc < length(testsets)) {
        if (showCompletion) {
            if (lc %in% printoutat) {
                cat(paste("\n *** ", round(100 * lc/length(testsets)), 
                  "% complete: tested ", lc, " of ", length(testsets), 
                  " sets of variables ", sep = ""))
            }
        }
        lc <- lc + 1
        usevariab <- testsets[[lc]]
        notusevariab <- (1:p)[-testsets[[lc]]]
        tmp <- getpval(Y, X[, c(1, 1 + usevariab), drop = FALSE], 
            ExpInd, test = test, maxNoObs = maxNoObs)
        pval <- tmp$pval
        Pall <- c(Pall, pval)
        if (pval > alpha) {
            if(stopIfEmpty){
                if(is.null(intersection)){
                    intersection <- usevariab
                }else{
                    intersection <- intersect(intersection,usevariab)
                }
                if( length(intersect)==0) cont <- FALSE
            }
                
            if (showAcceptedSets) 
                cat(paste("\n accepted set of variables ", paste(usevariab, 
                  collapse = ","), sep = ""))
            ConfInt[1, usevariab] <- pmax(ConfInt[1, usevariab, 
                drop = FALSE], tmp$coefficients + qnorm(1 - alpha/4) * 
                tmp$coefficientsvar, na.rm = TRUE)
            ConfInt[2, usevariab] <- pmin(ConfInt[2, usevariab, 
                drop = FALSE], tmp$coefficients - qnorm(1 - alpha/4) * 
                tmp$coefficientsvar, na.rm = TRUE)
            for (kc in usevariab) {
                Coeff[[kc]] <- c(Coeff[[kc]], tmp$coefficients[which(usevariab == 
                  kc)])
                CoeffVar[[kc]] <- c(CoeffVar[[kc]], tmp$coefficientsvar[which(usevariab == 
                  kc)])
            }
            if (length(notusevariab) >= 1) {
                ConfInt[1, notusevariab] <- pmax(ConfInt[1, notusevariab, 
                  drop = FALSE], 0, na.rm = TRUE)
                ConfInt[2, notusevariab] <- pmin(ConfInt[2, notusevariab, 
                  drop = FALSE], 0, na.rm = TRUE)
            }
            pvall <- c(pvall, pval)
        }
    }
    colnames(ConfInt) <- colnames(X[, -1, drop = FALSE])
    if (is.null(colnames(ConfInt))) 
        colnames(ConfInt) <- paste("Variable", 1:ncol(ConfInt))
    sig <- apply(sign(ConfInt[2:1, , drop = FALSE]), 2, function(x) prod(x))
    sigo <- sign(ConfInt[1, ])
    maximin <- sigo * apply(abs(ConfInt[2:1, , drop = FALSE]), 
        2, min) * (sig >= 0)
    pvalues <- rep(1, p)
    modelReject <- (!any(c(pvalempty,Pall) > gof))
    if (!modelReject) {
        for (k in usedvariables) {
            sel <- which(sapply(testsets, function(x, z) z %in% 
                x, k))
            pvalues[k] <- max(if (length(sel) > 0) max(Pall[-sel], 
                pvalempty) else max(Pall, pvalempty))
        }
    }else{
        ConfInt <-  NULL
        maximin <- NULL
    }
        
    retobj <- list(ConfInt = ConfInt, maximinCoefficients = maximin, 
        alpha = alpha, colnames = colnames(ConfInt), factor = is.factor(Y), 
        dimX = dim(X), Coeff = Coeff, CoeffVar = CoeffVar, modelReject = modelReject, 
        usedvariables = usedvariables, pvalues = pvalues, stopIfEmpty=stopIfEmpty, noEnv = length(ExpInd), gof=gof, bestModel=max(c(pvalempty,Pall)))
    class(retobj) <- "InvariantCausalPrediction"
    return(retobj)
}

