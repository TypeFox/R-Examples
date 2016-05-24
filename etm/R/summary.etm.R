summary.etm <- function(object, all = FALSE, ci.fun = "linear", level = 0.95, ...) {
    if (!inherits(object, "etm"))
        stop("'object' must be of class 'etm'")
    if (is.null(object$time)) {
        res <- list(P = object$est, s = object$s, t = object$t)
        class(res) <- "summary.etm"
        return(res)
    }
    if (level <= 0 | level > 1) {
        stop ("'level' must be between 0 and 1")
    }
    ref <- c("linear", "log", "cloglog", "log-log")
    if (sum(ci.fun %in% ref == FALSE) != 0) {
        stop("'ci.fun' is not correct. See help page")
    }
    if (all) {
        ind <- object$est != 0
        indi <- apply(ind, c(1, 2), function(temp){all(temp == FALSE)})
        tmp <- which(indi == FALSE, arr.ind = TRUE)
        tmp <- tmp[order(tmp[, 1]), ]
        namen <- list(rownames(indi), colnames(indi))
        trs <- lapply(seq_len(NROW(tmp)), function(i) {
            paste(namen[[1]][tmp[i, 1]], namen[[2]][tmp[i, 2]], sep = " ")
        })
        trs <- cbind(trs)
        absorb <- setdiff(levels(object$tran$to), levels(object$trans$from))
        for (i in seq_along(absorb))
            trs <- trs[-grep(paste("^", absorb[i], sep =""), trs, perl = TRUE)]
    }
    else {
        dtrs <- diag(outer(object$state.names, object$state.names, paste))
        absorb <- setdiff(levels(object$tran$to), levels(object$trans$from))
        for (i in seq_along(absorb))
            dtrs <- dtrs[-grep(paste("^", absorb[i], sep =""), dtrs, perl = TRUE)]
        tmp <- paste(object$trans[, 1], object$trans[, 2])
        trs <- c(tmp, dtrs)
    }
    res <- ci.transfo(object, trs, level, ci.fun)
    class(res) <- "summary.etm"
    res
}
        
