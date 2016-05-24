### plots.R:  Plot functions
### $Id: plots.R 42 2009-11-15 16:23:01Z bhm $

###
### Plot method for lspls objects
###

plot.lspls <- function(x, plottype = c("scores", "loadings"), ...) {
    plottype <- match.arg(plottype)
    plotFunc <- switch(plottype,
                       scores = scoreplot.lspls,
                       loadings = loadingplot.lspls)
    plotFunc(x, ...)
}


###
### Scoreplot
###

scoreplot.lspls <- function(object, ...) {
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(ask = TRUE)
    matnames <- strsplit(attr(object$terms, "term.labels")[-1], ":")
    for (i in seq(along = object$scores)) {
        if (is.matrix(object$scores[[i]])) {
            scoreplot(object$scores[[i]], comps = 1:object$ncomp[[i]],
                      main = matnames[[i]][1], ...)
        } else {
            for (j in seq(along = object$scores[[i]])) {
                scoreplot(object$scores[[i]][[j]],
                          comps = 1:object$ncomp[[i]][j],
                          main = matnames[[i]][j], ...)
            }
        }
    }
}


###
### Loadingplot
###

loadingplot.lspls <- function(object, ...) {
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mfrow = n2mfrow(length(unlist(object$ncomp))))
    matnames <- strsplit(attr(object$terms, "term.labels")[-1], ":")
    for (i in seq(along = object$loadings)) {
        if (is.matrix(object$loadings[[i]])) {
            loadingplot(object$loadings[[i]], comps = 1:object$ncomp[[i]],
                        main = matnames[[i]][1], ...)
        } else {
            for (j in seq(along = object$loadings[[i]])) {
                loadingplot(object$loadings[[i]][[j]],
                            comps = 1:object$ncomp[[i]][j],
                            main = matnames[[i]][j], ...)
            }
        }
    }
}


###
### Plot method for lsplsCv objects:
###

plot.lsplsCv <- function(x, which = c("RMSEP", "MSEP", "R2"), ncomp,
                         separate = TRUE, scale = !isTRUE(separate), ...) {
    which <- match.arg(which)
    val <- do.call(which, list(object = x, scale = scale))
    if (!isTRUE(separate)) {
        ## Aggregate over the responses, but keep a dummy dimension for
        ## the response (it simplifies the code below):
        dims <- c(1, dim(val)[-1])
        dns <- c(resp = "all responses", dimnames(val)[-1])
        if (which == "R2") {
            val <- array(colMeans(val), dim = dims, dimnames = dns)
        } else if (which == "RMSEP") {
            ## Dirty hack to get sqrt(sum(MSEP)) instead of sum(sqrt(MSEP)):
            val <- do.call("MSEP", list(object = x, scale = scale))
            val <- array(sqrt(colSums(val)), dim = dims, dimnames = dns)
        } else {
            val <- array(colSums(val), dim = dims, dimnames = dns)
        }
    }
    comps <- expand.grid(lapply(dimnames(val)[-1], as.numeric))
    ncomps <- rowSums(comps)
    ncombs <- nrow(comps)
    complabels <- apply(comps, 1, paste, collapse = "")
    mXlab <- if(missing(ncomp)) "total number of components" else "matrix"
    mYlab <- if (isTRUE(scale)) paste(which, "(std. resp.)") else which
    nResp <- dim(val)[1]
    if (nResp > 1) {
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = n2mfrow(nResp), oma = c(1, 1, 0, 0) + 0.1,
            mar = c(3, 3, 3, 1) + 0.1)
        xlab <- ""
        ylab <- ""
    } else {
        xlab <- mXlab
        ylab <- mYlab
    }
    respnames <-  dimnames(val)[[1]]
    if (missing(ncomp)) {
        val <- aperm(val, c(2:length(dim(val)), 1)) # Make "resp" the last dimension
        for (i in 1:nResp) {
            cval <- c(val)[ncombs * (i - 1) + 1:ncombs]
            plot(ncomps, cval, type = "n", xlab = xlab, ylab = ylab,
                 main = respnames[i], ...)
            text(ncomps, cval, labels = complabels)
            oncomps <- min(ncomps):max(ncomps)
            bestval <- numeric(length(oncomps))
            for (i in seq(along = oncomps))
                bestval[i] <- if (which == "R2") max(cval[ncomps == oncomps[i]])
                              else min(cval[ncomps == oncomps[i]])
            lines(oncomps, bestval, lty = 2, col = 2)
        } ## for
    } else {
        ## Extract a matrix with measure values versus the matrices included,
        ## for the specified number of components
        nMat <- length(ncomp) + 1
        matNames <- attr(terms(x), "term.labels")
        if(nMat != length(matNames))
            stop("'ncomp' must contain ", length(matNames) - 1, " elements")
        plotVals <- matrix(nrow = nResp, ncol = nMat)
        valInds <- as.list(rep(1, length(unlist(ncomp)))) # Indices into val
        plotVals[,1] <- do.call("[", c(list(val, 1:nResp), valInds))
        l <- 0                          # index into valInds
        for (j in seq_along(ncomp)) {
            for(k in seq_along(ncomp[[j]])) {
                l <- l + 1
                valInds[[l]] <- ncomp[[j]][k] + 1
            }
            plotVals[,j+1] <- do.call("[", c(list(val, 1:nResp), valInds))
        }
        for (i in 1:nResp) {
            plot(plotVals[i,], type = "b", xlab = xlab, ylab = ylab,
                 main = respnames[i], xaxt = "n", ...)
            axis(1, at = seq_along(matNames), labels = matNames)
        }
    } # if (missing(ncomp)) ... else
    if (nResp > 1) {
        ## Add outer margin text:
        mtext(mXlab, side = 1, outer = TRUE)
        mtext(mYlab, side = 2, outer = TRUE)
    }
} ## function
