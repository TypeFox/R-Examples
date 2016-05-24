`boxplot.mefa` <-
function(x, stat=1:4, all = TRUE, show=TRUE, ylab, xlab, ...)
{
    if (is.null(x$segm) || dim(x)[3] == 1)
        stop("at least 2 segments needed")
    if (!all(stat %in% 1:4))
        stop("'stat' must be in 1:4")
    if (!length(stat) == 1) stat <- 1
    if (all)
        k <- 1 else k <- 0
    yval <- list()
    if (all)
        yval[[1]] <- summary(x)[[stat]]
    for (i in (1 + k):(dim(x)[3] + k))
        yval[[i]] <- summary(mefa(x$segm[[(i - k)]]))[[stat]]
    yval <- unlist(yval)
    if (all) {
    levs <- if ("all" %in% dimnames(x)$segm)
        paste("segm", dimnames(x)$segm, sep=".") else dimnames(x)$segm
    levs <- tolower(c("all", levs))
    } else levs <- tolower(dimnames(x)$segm)
    xval <- as.factor(rep(levs, each=length(summary(x)[[stat]])))
    if (missing(ylab))
        ylab <- c("Frequency of taxa", "Frequency of individuals",
        "Frequency of occurrence", "Abundance")[stat]
    if (missing(xlab))
        xlab <- "Segments"
    if (show)
        boxplot(yval ~ xval, xlab=xlab, ylab=ylab, ...)
    if (show)
# invisibly returns plotted walues
        invisible(cbind(x=xval, y=yval)) else return(cbind(x=xval, y=yval))
}

