reduc <-
function (x, clu, labels = NULL) 
{
    if (isTRUE(is.array(x) == TRUE) == FALSE) 
        stop("'x' must be an array object.")
    if (isTRUE(length(clu) != dim(x)[1]) == TRUE) 
        stop("'clu' does not match the order of 'x'.")
    if (is.character(clu) == FALSE) {
        ifelse(isTRUE(0L %in% as.numeric(levels(factor(clu)))) == 
            TRUE, clu <- clu + 1L, NA)
    }
    else if (is.character(clu) == TRUE) {
        tmp <- clu
        for (i in 1:nlevels(factor(clu))) {
            clu[which(levels(factor(tmp))[i] == clu)] <- i
        }
        rm(i)
        clu <- methods::as(clu, "numeric")
        rm(tmp)
    }
    else {
        NA
    }
    clu[which(is.na(clu))] <- max(as.numeric(levels(factor(clu)))) + 
        1L
    lngt <- nlevels(factor(clu))
    or <- list()
    for (i in as.numeric(levels(factor(clu)))) {
        or[[i]] <- which(clu == i)
    }
    rm(i)
    cls <- list()
    length(cls) <- lngt
    k <- 1L
    for (i in 1:length(or)) {
        if (isTRUE(is.null(or[[i]])) == FALSE) {
            cls[[k]] <- or[[i]]
            k <- k + 1L
        }
        else {
            NA
        }
    }
    rm(i, k)
    if (isTRUE(is.na(dim(x)[3]) == TRUE)) {
        if (isTRUE(is.null(dimnames(x)[[1]]) == TRUE) == TRUE) 
            dimnames(x)[[1]] <- dimnames(x)[[2]] <- 1:nrow(x)
        bm <- array(dim = c(lngt, lngt))
        for (i in 1:lngt) {
            for (j in 1:lngt) {
                bm[i, j] <- sum(x[cls[[i]], cls[[j]]])
            }
        }
        rm(i, j)
        bm <- dichot(bm)
        ifelse(is.null(labels) == FALSE, rownames(bm) <- colnames(bm) <- labels, 
            NA)
        return(bm)
    }
    else if (isTRUE(is.na(dim(x)[3]) == FALSE)) {
        px <- x
        bm <- array(dim = c(lngt, lngt, dim(x)[3]))
        for (k in 1:dim(x)[3]) {
            for (i in 1:lngt) {
                for (j in 1:lngt) {
                  bm[i, j, k] <- sum(px[cls[[i]], cls[[j]], k])
                }
            }
            rm(i, j)
        }
        rm(k)
        bm <- dichot(bm, c = 1L)
        if (is.null(labels) == FALSE) 
            dimnames(bm)[[1]] <- dimnames(bm)[[2]] <- labels
        if (is.null(dimnames(x)[[3]]) == FALSE) 
            dimnames(bm)[[3]] <- dimnames(x)[[3]]
        return(bm)
    }
}
