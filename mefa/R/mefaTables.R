`mefaTables` <-
function(xtab, dframe, margin, index=NULL, drop.index=FALSE, xtab.fixed=TRUE)
{
    if (margin != 1 && margin != 2)
        stop("'margin' should be 1 or 2")
    if (NCOL(dframe) == 1) {
        dframe <- data.frame(x=dframe, copy=dframe)
        onecol <- TRUE
    } else onecol <- FALSE
# ordering tables
    if (margin == 1) {
        rank.orig <- c(1:nrow(xtab))[order(rownames(xtab))]
        nam.orig <- rownames(xtab)
        xtab <- xtab[order(rownames(xtab)), ]
        xnam <- rownames(xtab)
        } else {
        rank.orig <- c(1:ncol(xtab))[order(colnames(xtab))]
        nam.orig <- colnames(xtab)
        xtab <- xtab[, order(colnames(xtab))]
        xnam <- colnames(xtab)}
# if index is specified that used as rownames
    if (!is.null(index)) 
        rownames(dframe) <- dframe[, index]
    dnam <- rownames(dframe)
# xtab fixed
    if (xtab.fixed) {
        dsub <- dframe[dnam %in% xnam, ]
        dsub <- dsub[order(rownames(dsub)), ]
        xsub <- xtab
        rank.final <- rank.orig
# xtab not fixed
    } else {
        int <- intersect(dnam, xnam)
        dsub <- dframe[dnam %in% int, ]
        dsub <- dsub[order(rownames(dsub)), ]
        rank.final <- rank.orig[nam.orig %in% int]
        xsub <- if (margin == 1)
            xtab[xnam %in% int, ] else xsub <- xtab[, xnam %in% int]
        }
# here ordering again to get back original
    xsub <- if (margin == 1)
        xsub[order(rank.final), ] else xsub[, order(rank.final)]
    dsub <- dsub[order(rank.final), ]
    nsub <- if (margin == 1)
        rownames(xsub) else colnames(xsub)
# check result
    if (!identical(rownames(dsub), nsub))
        stop("names do not match")
# drop index
    if (!is.null(index) && drop.index)
        dsub[, index] <- NULL
    if (onecol)
        dsub <- as.data.frame(x=dsub[,1])
    return(list(xtab=xsub, dtab=dsub))
}

