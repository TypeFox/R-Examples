cph <-
function (W, labels = NULL) 
{
    if (isTRUE(attr(W, "class") == "Rel.Box") == FALSE) 
        stop("\"W$W\" must be a \"Rel.Box\" class.")
    rls <- array(dim = c(dim(W$W)[3], dim(W$W)[1], dim(W$W)[1]))
    for (i in 1:dim(W$W)[3]) {
        for (j in 1:dim(W$W)[1]) {
            rls[i, , j] <- W$W[j, , i]
        }
    }
    phs <- array(0, dim = c(dim(W$W)[1], dim(W$W)[1], dim(W$W)[1]))
    for (k in 1:dim(W$W)[1]) {
        for (j in 1:dim(W$W)[1]) {
            for (i in 1:dim(W$W)[1]) {
                if ((as.numeric(any(rls[, i, k] < rls[, j, k])) == 
                  1 && as.numeric(any(rls[, j, k] < rls[, i, 
                  k])) == 0) | as.numeric(all(rls[, i, k] == 
                  rls[, j, k])) == 1) 
                  phs[i, j, k] <- 1
            }
        }
    }
    for (k in 1:dim(W$W)[1]) {
        for (i in 1:dim(W$W)[1]) {
            if (sum(rls[, i, k]) == 0) 
                phs[i, , k] <- 0
        }
    }
    dimnames(phs)[[1]] <- dimnames(phs)[[2]] <- dimnames(W$W)[[1]]
    tmp <- data.frame(matrix(ncol = (dim(phs)[1] * dim(phs)[2]), 
        nrow = 0))
    for (i in 1:dim(phs)[3]) {
        ifelse(dim(phs)[3] > 1, tmp[i, ] <- as.vector(phs[, , 
            i]), tmp <- as.vector(phs))
    }
    rm(i)
    phu <- matrix(replace(sapply(tmp, sum), sapply(tmp, sum) >= 
        1, 1), nrow = dim(W$W)[1], ncol = dim(W$W)[1])
    diag(phu) <- 1
    for (i in seq_len(ncol(phu))) {
        tmp <- outer(phu[, i], phu[i, ], pmin.int)
        phu <- pmax(phu, tmp)
    }
    if (isTRUE(is.null(labels)) == FALSE) {
        ifelse(isTRUE(length(labels) == dim(phu)[1]) == TRUE, 
            dimnames(phu)[[1]] <- dimnames(phu)[[2]] <- labels, 
            NA)
    }
    else if (isTRUE(is.null(labels)) == TRUE) {
        dimnames(phu)[[1]] <- dimnames(phu)[[2]] <- W$lbs
    }
    return(phu)
}
