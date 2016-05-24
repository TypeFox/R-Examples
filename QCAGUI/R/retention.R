`retention` <- 
function(data, outcome = "", conditions = "", type = "corruption", 
         dep = TRUE, n.cut = 1, incl.cut = 1, p.pert = 0.5, n.pert = 1) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    names(data) <- toupper(names(data))
    conditions <- toupper(conditions)
    outcome <- toupper(outcome)
    
    if (identical(conditions, "")) {
        conditions <- names(data)[-which(names(data) == outcome)]
    }
    else {
        conditions <- QCA::splitstr(conditions)
    }
    
    data <- data[, c(conditions, outcome)]
    udata <- unique(data[, conditions])
    rownames(udata) <- seq(nrow(udata))
    cpos <- cneg <- rep(0, nrow(udata))
    for (i in seq(nrow(udata))) {
        for (j in seq(nrow(data))) {
            if (all(udata[i, ] == data[j, conditions])) {
                if (data[j, outcome] == 1) {
                    cpos[i] <- cpos[i] + 1
                }
                else if (data[j, outcome] == 0) {
                    cneg[i] <- cneg[i] + 1
                }
            }
        }
    }
    
    total <- cpos + cneg
    udata <- udata[total >= n.cut, , drop = FALSE]
    cpos <- cpos[total >= n.cut]
    cneg <- cneg[total >= n.cut]
    total <- total[total >= n.cut]
    
    calculatePairs <- function(x, n.pert, type = "deletion") {
        pairsxl <- combn(nrow(udata), min(x, nrow(udata)))
        nofsetsxl <- 0
        for (j in seq(ncol(pairsxl))) {
            cposneg <- NULL
            for (l in seq(min(x, nrow(pairsxl)))) {
                cposneg <- c(cposneg, cbind(cpos, cneg)[pairsxl[l, j], ])
            }
            allpairs <- createMatrix(cposneg + 2) - 1
            allpairs <- allpairs[apply(allpairs, 1, function(y) all(y >= 0)), , drop = FALSE]
            linesubset <- rep(TRUE, nrow(allpairs))
            for (l in seq(min(x, nrow(pairsxl)))) {
                linesubset <- linesubset & rowSums(allpairs[, l*2 - c(1, 0)]) >= 1
            }
            allpairs <- allpairs[linesubset & rowSums(allpairs) <= n.pert, , drop = FALSE]
            for (i in seq(nrow(allpairs))) {
                lchanges <- rep(FALSE, min(x, nrow(pairsxl)))
                for (l in seq(min(x, nrow(pairsxl)))) {
                    initially <- cpos[pairsxl[l, j]]/total[pairsxl[l, j]]
                    if (type == "deletion") {
                        newtotaless <- total[pairsxl[l, j]] - allpairs[i, l*2 - 1]
                        after <- (cpos[pairsxl[l, j]] - allpairs[i, l*2 - 1])/newtotaless
                        lchanges[l] <- ((initially >= incl.cut & after <  incl.cut) | newtotaless <  n.cut) |
                                       ((initially <  incl.cut & after >= incl.cut) & newtotaless >= n.cut)
                    }
                    else if (type == "corruption") {
                        after <- (cpos[pairsxl[l, j]] + allpairs[i, l*2] - allpairs[i, l*2 - 1])/total[pairsxl[l, j]]
                        lchanges[l] <- (initially >= incl.cut & after <  incl.cut) |
                                       (initially <  incl.cut & after >= incl.cut)
                    }
                }
                if (all(lchanges)) {
                    combs <- 1
                    for (l in seq(min(x, nrow(pairsxl)))) {
                        combs <- combs*choose(cpos[pairsxl[l, j]], allpairs[i, l*2 - 1])
                        combs <- combs*choose(cneg[pairsxl[l, j]], allpairs[i, l*2])
                    }
                    nofsetsxl <- nofsetsxl + combs*choose(nrow(data) - sum(cposneg), n.pert - sum(allpairs[i, ]))
                }
            }
        }
        return(nofsetsxl)
    }
    
    
    if (dep) {
        nofsets <- 0
        totalsets <- choose(nrow(data), n.pert)
        for (i in seq(n.pert)) {
            if (nofsets != totalsets) {
                nofsetsxl <- calculatePairs(i, n.pert, type = type)
                nofsets <- nofsets + ifelse(i %% 2 == 1, nofsetsxl, -1*nofsetsxl)
            }
        }
        
        return(as.vector(1 - nofsets/totalsets))
    }
    else {
        pfinal <- 1
        if (type == "deletion") {
            for (l in seq(nrow(udata))) {
                ptmp <- 1
                allpairs <- createMatrix(c(cpos[l], cneg[l]) + 2) - 1
                allpairs <- allpairs[apply(allpairs, 1, function(x) all(x >= 0)), , drop = FALSE]
                allpairs <- allpairs[rowSums(allpairs) >= 1, , drop = FALSE]
                for (i in seq(nrow(allpairs))) {
                    newtotaless <- total[l] - allpairs[i, 1] - allpairs[i, 2]
                    initially <- cpos[l]/total[l]
                    after <- (cpos[l] - allpairs[i, 1])/newtotaless
                    if (((initially >= incl.cut & after <  incl.cut) | newtotaless <  n.cut) |
                        ((initially <  incl.cut & after >= incl.cut) & newtotaless >= n.cut)) {
                           ptmp <- ptmp - dbinom(allpairs[i, 1], cpos[l], p.pert) * dbinom(allpairs[i, 2], cneg[l], p.pert)
                    }
                }
                pfinal <- pfinal*ptmp
            }
        }
        else if (type == "corruption") {
            for (l in seq(nrow(udata))) {
                ptmp <- 1
                allpairs <- createMatrix(c(cpos[l], cneg[l]) + 2) - 1
                allpairs <- allpairs[apply(allpairs, 1, function(x) all(x >= 0)), , drop = FALSE]
                allpairs <- allpairs[rowSums(allpairs) >= 1, , drop = FALSE]
                for (i in seq(nrow(allpairs))) {
                    initially <- cpos[l]/total[l]
                    after <- (cpos[l] - allpairs[i, 1] + allpairs[i, 2])/total[l]
                    if ((initially >= incl.cut & after < incl.cut) | (initially < incl.cut & after >= incl.cut)) {
                           ptmp <- ptmp - dbinom(allpairs[i, 1], cpos[l], p.pert) * dbinom(allpairs[i, 2], cneg[l], p.pert)
                    }
                }
                pfinal <- pfinal*ptmp
            }
        }
        
        return(pfinal)
    }
}

