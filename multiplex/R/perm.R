perm <-
function (x, clu, rev = FALSE) 
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
    or <- list()
    for (i in as.numeric(levels(factor(clu)))) {
        or[[i]] <- which(clu == i)
    }
    rm(i)
    if (isTRUE(any(is.na(unlist(or)))) == TRUE) {
        for (i in 1:length(or)) {
            ifelse(isTRUE(is.null(or[[i]])) == TRUE, or[[i]] <- NA, 
                NA)
        }
        rm(i)
        nor <- as.vector(stats::na.omit(unlist(or)))
    }
    else {
        nor <- order(unlist(or))
    }
    if (isTRUE(is.na(dim(x)[3]) == TRUE)) {
        prm <- vector()
        for (i in nor) {
            prm[nor[i]] <- i
        }
        rm(i)
        ifelse(isTRUE(rev == TRUE) == TRUE, return(x[rev(prm), 
            rev(prm)]), return(x[prm, prm]))
    }
    else {
        px <- x
        for (k in 1:dim(px)[3]) {
            prm <- vector()
            for (i in nor) {
                prm[nor[i]] <- i
            }
            rm(i)
            px[, , k] <- px[prm, prm, k]
        }
        rm(k)
        if (isTRUE(is.null(dimnames(x)[[1]]) == FALSE)) {
            lbs <- vector()
            length(lbs) <- length(clu)
            for (i in 1:length(nor)) lbs[i] <- dimnames(x)[[1]][which(nor == 
                i)]
            dimnames(px)[[1]] <- dimnames(px)[[2]] <- lbs
        }
        ifelse(isTRUE(is.null(dimnames(x)[[3]]) == FALSE), dimnames(px)[[3]] <- dimnames(x)[[3]], 
            NA)
        return(px)
    }
}
