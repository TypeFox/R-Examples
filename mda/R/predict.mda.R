predict.mda <-
function (object, newdata, type = c("class", "variates", "posterior", 
    "hierarchical", "weights"), prior = NULL, dimension = R - 
    1, g, ...) 
{
    type <- match.arg(type)
    Rk <- dim(object$means)
    R <- Rk[1]
    k <- Rk[2]
    if (type == "hierarchical") {
        if (missing(dimension)) 
            dimension.set <- seq(k)
        else {
            dimension.set <- dimension[dimension <= k]
            if (!length(dimension.set)) 
                dimension.set <- k
            dimension <- max(dimension.set)
        }
        Pclass <- vector("list", length(dimension.set))
        names(Pclass) <- paste("D", dimension.set, sep = "")
        for (ad in seq(along = dimension.set)) {
            d <- dimension.set[ad]
            Pclass[[ad]] <- if (missing(newdata)) 
                Recall(object, prior = prior, dimension = d, 
                  ...)
            else Recall(object, newdata, prior = prior, dimension = d, 
                ...)
        }
        rownames <- names(Pclass[[1]])
        if (is.null(rownames)) 
            rownames <- paste(seq(along = Pclass[[1]]))
        return(structure(Pclass, class = "data.frame", row.names = rownames, 
            dimensions = dimension.set))
    }
    else dimension <- min(max(dimension), k)
    if (is.null(prior)) 
        prior <- object$prior
    else {
        if (any(prior < 0) | round(sum(prior), 5) != 1) 
            stop("innappropriate prior")
    }
    if (type == "variates") 
        return(NextMethod("predict"))
    rowmin <- function(mat) {
        ncc <- ncol(mat)
        if (ncc == 1) 
            return(drop(mat))
        rowm <- pmin(mat[, 1], mat[, 2])
        if (ncc == 2) 
            return(rowm)
        else {
            for (j in seq(from = 3, to = ncc)) rowm <- pmin(rowm, 
                mat[, j])
        }
        rowm
    }
    dmat <- if (missing(newdata)) 
        predict.fda(object, type = "distances", dimension = dimension, 
            ...)
    else predict.fda(object, newdata, type = "distances", dimension = dimension, 
        ...)
    Assign <- object$assign
    sub.prior <- object$sub.prior
    J <- length(Assign)
    if (type == "weights") {
        if (missing(newdata)) 
            return(object$weights)
        g <- as.numeric(g)
        weights <- Assign
        for (j in seq(J)) {
            TT <- dmat[g == j, Assign[[j]], drop = FALSE]
            TT <- exp(-0.5 * (TT - rowmin(TT)))
            TT <- TT * outer(rep(1, nrow(TT)), sub.prior[[j]])
            weights[[j]] <- TT/drop(TT %*% rep(1, ncol(TT)))
        }
        return(weights)
    }
    pclass <- matrix(1, nrow(dmat), J)
    dmat <- exp(-0.5 * (dmat - rowmin(dmat)))
    for (j in seq(J)) {
        TT <- dmat[, Assign[[j]], drop = FALSE]
        TT <- TT * outer(rep(1, nrow(TT)), sub.prior[[j]])
        pclass[, j] <- prior[j] * drop(TT %*% rep(1, ncol(TT)))
    }
    dimnames(pclass) <- list(NULL, names(Assign))
    switch(type, class = softmax(pclass), posterior = pclass/drop(pclass %*% 
        rep(1, J)))
}

