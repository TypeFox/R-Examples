`join` <- function(..., verbose = FALSE, na.replace = TRUE,
                   split = TRUE, value = 0,
                   type = c("outer", "left", "inner")) {
    outerJoin <- function(X) {
        ## From code provided by Sundar Dorai-Raj in R-Help posting:
        ## http://article.gmane.org/gmane.comp.lang.r.general/63042/
        cn <- unique(unlist(lapply(X, colnames)))
        for(i in seq(along = X)) {
            if(any(m <- !cn %in% colnames(X[[i]]))) {
                na <- matrix(NA, nrow(X[[i]]), sum(m))
                dimnames(na) <- list(rownames(X[[i]]), cn[m])
                X[[i]] <- cbind(X[[i]], na)
            }
        }
        joined <- do.call("rbind", X)
        colnames(joined) <- cn
        joined
    }
    leftJoin <- function(X) {
        cn <- unique(unlist(lapply(X, colnames)[[1]]))
        ## if more than 2 df in X, merge all bar first
        if(length(X) > 2)
            dfs <- outerJoin(X[-1])
        else
            dfs <- X[[2]]
        ## matched column names
        mcn <- match(colnames(dfs), cn)
        mcn2 <- match(cn, colnames(dfs))
        mcn <- mcn[!is.na(mcn)]
        mcn2 <- mcn2[!is.na(mcn2)]
        joined <- matrix(NA, ncol = dims[1,2], nrow = sum(dims[,1]))
        joined[1:dims[1,1], ] <- data.matrix(X[[1]])
        joined[(dims[1,1]+1):NROW(joined), mcn] <- data.matrix(dfs[, mcn2])
        colnames(joined) <- cn
        joined
    }
    innerJoin <- function(X) {
        cn <- lapply(X, colnames)
        cn <- Reduce(intersect, cn)
        ##joined <- matrix(NA, ncol = length(cn), nrow = sum(dims[,1]))
        joined <- vector(length = length(X), mode = "list")
        for(i in seq_along(joined)) {
            joined[[i]] <- data.matrix(X[[i]][, cn])
        }
        joined <- do.call(rbind, joined)
        colnames(joined) <- cn
        joined
    }
    x <- list(...)
    if(any(!sapply(x, inherits, "data.frame", USE.NAMES = FALSE)))
        stop("\nall objects to be merged must be data frames.")
    dims <- t(sapply(x, dim))
    n.joined <- nrow(dims)
    if(missing(type))
        type <- "outer"
    type <- match.arg(type)
    joined <- switch(type,
                     outer = outerJoin(x),
                     left = leftJoin(x),
                     inner = innerJoin(x))
    if(na.replace) {
        joined[is.na(joined)] <- value
    }
    rn <- lapply(x, rownames)
    if(verbose) {
        stats <- rbind(dims, dim(joined))
        rownames(stats) <- c(paste("Data set ", c(1:n.joined), ":", sep = ""),
                             "Merged:")
        colnames(stats) <- c("Rows", "Cols")
        cat("\nSummary:\n\n")
        printCoefmat(stats, digits = max(3, getOption("digits") - 3),
                     na.print = "")
        cat("\n")
    }
    if(split) {
        retval <- vector(mode = "list", length = n.joined)
        ends<- cumsum(dims[,1])
        start <- c(1, ends[-n.joined] + 1)
        for(i in seq_len(n.joined)) {
            retval[[i]] <- as.data.frame(joined[start[i]:ends[i], ])
            rownames(retval[[i]]) <- rn[[i]]
        }
        names(retval) <- as.character(match.call())[2:(n.joined+1)]
        class(retval) <- "join"
    } else {
        retval <- as.data.frame(joined, row.names = rownames(joined))
        class(retval) <- c("join", class(retval))
    }
    attr(retval, "type") <- type
    retval
}
