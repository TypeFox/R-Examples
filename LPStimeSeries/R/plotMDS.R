plotMDS <- function(object, newdata, classinfo, k=2, palette=NULL, pch=20, ...) {
    if (!inherits(object, "learnPattern")) 
        stop(deparse(substitute(object)), " must be a learnPattern object")

	if(!is.factor(classinfo)) classinfo <- as.factor(classinfo)
    op <- par(pty="s")
    on.exit(par(op))

    sim <- computeSimilarity(object,newdata,newdata)
    sim <- sim/(2*sum(object$nobs))	
    lps.mds <- stats::cmdscale(sim, eig=TRUE, k=k)
    colnames(lps.mds$points) <- paste("Dim", 1:k)
    nlevs <- nlevels(classinfo)
    if (is.null(palette)) {
        palette <- if (nlevs < 12)
            RColorBrewer::brewer.pal(nlevs, "Set1") else rainbow(nlevs)
    }
    if (k <= 2) {
        plot(lps.mds$points, col=palette[as.numeric(classinfo)], pch=pch, ...)
    } else {
        pairs(lps.mds$points, col=palette[as.numeric(classinfo)], pch=pch, ...)
    }
    invisible(lps.mds)
}
