`[.mefa` <-
function (x, i=1:dim(x)[1], j=1:dim(x)[2], k=1:dim(x)[3], drop=FALSE)
{
# check i,j,k values
    if (length(i) == 0 || length(j) == 0 || length(k) == 0)
        return(x)

    if (is.character(i))
        i <- match(i, dimnames(x)[[1]])
    if (is.character(j))
        j <- match(j, dimnames(x)[[2]])
    if (is.character(k))
        k <- match(k, dimnames(x)[[3]])

    if ((any(i < 0) & !all(i < 0)) || (any(j < 0) & !all(j < 0)) || (any(k < 0) & !all(k < 0)))
        stop("unresolvable indexing")
    if (all(i < 0))
        i <- which(!(1:dim(x)[1] %in% abs(i)))
    if (all(j < 0))
        j <- which(!(1:dim(x)[2] %in% abs(j)))
    if (all(k < 0))
        k <- which(!(1:dim(x)[3] %in% abs(k)))

    if (length(i) == 1 || length(j) == 1)
        stop("dimensions are too few to be extracted")

# logical if there is a need to extract segments based on k values
    subsegm <- if (!is.null(x$segm))
        length(k) != length(x$segm) else FALSE
    x$call <- match.call()
# here segments are extracted based on i and j
    if (is.null(x$segm)) {
        segm.out <- NULL} else {
        segm.out <- list()
        for (nn in 1:length(k))
            segm.out[[nn]] <- x$segm[[k[nn]]][i, j]
        names(segm.out) <- names(x$segm)[k]
        }
# xtab is deifferent when segments are nested
    if (!subsegm) {
        xtab.out <- x$xtab[i, j]
        } else {
        if (!is.null(x$segm) && length(k) == 1) {
            xtab.out <- segm.out[[1]]
            } else {
                if (attr(x, "nested")) {
                    xtab.out <- x$segm[[max(k)]]
                    } else {
                    xtab.out <- segm.out[[1]]
                    for (nn in 2:length(k))
                        xtab.out <- xtab.out + segm.out[[nn]] }
            }
        }
# correction of nested segment names
    if (subsegm && attr(x, "nested")) {
        mn <- match(names(x$segm)[k], names(x$segm))
        if (!identical(mn, mn[order(mn)]))
            stop("nested segments cannot be reordered")
        nam1 <- names(x$segm)[1]
        nam2 <- names(x$segm)[names(x$segm) %in% names(x$segm)[k]][1]
        newnam <- names(x$segm)[k]
        if (!identical(nam1, nam2)) {
            for (i in 1:(nchar(nam1)+1)) {
                substr(nam2, 1, 1) <- ""
                substr(newnam, 1, 1) <- ""
            }
        names(segm.out) <- if (length(k) == 1)
            nam2 else c(nam2, paste(nam2, newnam[-1], sep = "-"))
        }
    }
# put together the psrts
    x$segm <- segm.out
    x$xtab <- xtab.out
# extract samp and taxa tables and apply drop argument
    if (is.null(x$samp))
        x$samp <- NULL else {
        if (!is.data.frame(x$samp))
            stop("'x$samp' is not a data frame")
        x$samp <- x$samp[i, ]
        if (drop) x$samp[] <- lapply(x$samp, function(x) x[drop = TRUE])}
    if (is.null(x$taxa))
        x$taxa <- NULL else {
        if (!is.data.frame(x$taxa))
            stop("'x$samp' is not a data frame")
        x$taxa <- x$taxa[j, ]
        if (drop) x$taxa[] <- lapply(x$taxa, function(x) x[drop = TRUE])}
    if (subsegm && length(k) == 1)
        attr(x, "nested") <- FALSE
    return(x)
}

