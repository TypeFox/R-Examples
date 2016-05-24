`aggregate.mefa` <-
function(x, by.samp=NULL, by.taxa=NULL, ...)
{
    if (is.null(by.samp) && is.null(by.taxa))
        return(x) else {

    xtab <- as.data.frame(x$xtab)

# aggregation based on sample indices
    if (!is.null(by.samp)) {
        if (!is.object(by.samp))
            if (length(by.samp) == 1)
                by.samp <- x$samp[, by.samp] else {
                by.samp <- interaction(x$samp[, by.samp])}
        if (length(unique(by.samp)) == 1)
            stop("'by.samp' should contain at least 2 levels")
# samples table is NULL because the FUN of aggregation is not trivial
        x$samp <- NULL
        xtab <- aggregate(xtab, list(by.samp), sum, ...)
        rownames(xtab) <- xtab[,1]
        xtab[,1] <- NULL
# aggregation of segments
        if (!is.null(x$segm)){
            for (i in 1:length(x$segm)) {
                x$segm[[i]] <- aggregate(x$segm[[i]], list(by.samp), sum, ...)
                rownames(x$segm[[i]]) <- x$segm[[i]][,1]
                x$segm[[i]][,1] <- NULL}}
        }

# aggregation based on taxon indices
    if (!is.null(by.taxa)) {
        if (!is.object(by.taxa))
            if (length(by.taxa) == 1)
                by.taxa <- x$taxa[, by.taxa] else {
                by.taxa <- interaction(x$taxa[, by.taxa])}
        if (length(unique(by.taxa)) == 1)
            stop("'by.taxa' should contain at least 2 levels")
# taxa table is NULL because the FUN of aggregation is not trivial
        x$taxa <- NULL
        xtab <- aggregate(t(xtab), list(by.taxa), sum, ...)
        rownames(xtab) <- xtab[,1]
        xtab[,1] <- NULL
        xtab <- t(xtab)
# aggregation of segments
        if (!is.null(x$segm)){
            for (i in 1:length(x$segm)) {
                x$segm[[i]] <- aggregate(t(x$segm[[i]]), list(by.taxa), sum, ...)
                rownames(x$segm[[i]]) <- x$segm[[i]][,1]
                x$segm[[i]][,1] <- NULL
                x$segm[[i]] <- t(x$segm[[i]])}}
        }
    x$xtab <- as.matrix(xtab)
    x$call <- match.call()
    return(x)
    }
}

