
# $Id: maxsets.R 229 2008-04-04 11:47:30Z thothorn $

### compute all possible (ordered) subsets of the index set K
### cf. Westfall (1997, Section 3.2)
allsubsets <- function(K) 
{
    if (length(K) == 0) return(list(NULL))
    if (length(K) == 1) return(list(K))
    ret <- as.list(K)
    for (i in 1:(length(K)-1)) {
        tmp <- allsubsets(K[-(1:i)])
        for (j in 1:length(tmp))
            tmp[[j]] <- c(K[i], tmp[[j]])
        ret <- c(ret, tmp)
    }
    ret
}

### check if any of C[,1:(min(K)-1)] is in column space of C[,K]
### cf. Westfall (1997, Section 3.2)
checkCS <- function(K, C) 
{
    if (length(K) == ncol(C)) return(TRUE)
    CK <- C[,K,drop = FALSE]
    Cj <- C[,1:(min(K)-1), drop = FALSE]
    tmp <- Cj - (CK %*% MPinv(CK)$MPinv %*% Cj)
    all(colSums(tmp^2) > .Machine$double.eps)
}

### remove redundant index sets
rmsets <- function(sets) 
{
    if (length(sets) == 1) return(sets)
    rm <- logical(length(sets))
    for (j in 1:(length(sets) - 1)) {
        set <- sets[[j]]
        rm[j] <- any(sapply(sets[(j + 1):length(sets)], function(x) 
                            all(set %in% x)))
    }
    sets[!rm]
}

### compute maximal sets of linear hypotheses
### cf. Westfall (1997, Section 3.2)
maxsets <- function(K) 
{
    C <- t(K)
    k <- ncol(C)
    S <- 1:k
    ret <- vector(mode = "list", length = k)

    for (j in S) {
        tmp <- allsubsets(S[-(1:j)])
        for (i in 1:length(tmp))
            tmp[[i]] <- c(j, tmp[[i]])
        if (length(tmp) > 1 || length(tmp[[1]]) > 1)
            tmp <- c(j, tmp)
        ret[[j]] <- tmp[sapply(tmp, checkCS, C = C)]
        
    }
    lapply(ret, rmsets)
}
