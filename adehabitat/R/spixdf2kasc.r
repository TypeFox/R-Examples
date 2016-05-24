spixdf2kasc <- function (sg)
{
    if (inherits(sg, "SpatialPixelsDataFrame"))
        sg <- as(sg, "SpatialGridDataFrame")
    if (!inherits(sg, "SpatialGridDataFrame"))
        stop(paste("sg should be of class \"SpatialPixelsDataFrame\"",
            "\nor \"SpatialGridDataFrame\""))
    gr <- gridparameters(sg)
    if (nrow(gr) > 2)
        stop("sg should be defined in two dimensions")
    if (gr[1, 2] != gr[2, 2])
        stop("the cellsize should be the same in x and y directions")
    nRow <- gr[2, 3]
    nRows <- nRow:1
    nCol <- gr[1, 3]
    nCols <- 1:nCol
    lo <- integer(length=nRow*nCol)
    for (i in 1:nRow) {
        r <- ((i-1)*(nCol)+1):(i*nCol)
        lo[r] <- (nRows[i]-1)*(nCol) + nCols
    }
    uu <- names(sg)
    if (length(uu) ==1) {
        ka <- data.frame(slot(sg, "data")[lo,])
        names(ka) <- uu
    } else ka <- slot(sg, "data")[lo,]
    attr(ka, "cellsize") <- gr[2, 2]
    attr(ka, "xll") <- gr[1, 1]
    attr(ka, "yll") <- gr[2, 1]
    attr(ka, "ncol") <- gr[1, 3]
    attr(ka, "nrow") <- gr[2, 3]
    class(ka) <- c("kasc", "data.frame")
    if (ncol(ka) == 1) {
        v <- ka[, 1]
        if ((is.numeric(v)) | (is.logical(v))) {
            e <- matrix(v, ncol = attr(ka, "nrow"))
            attr(e, "type") <- "numeric"
        }
        else {
            tc2 <- levels(v)
            v <- as.numeric(v)
            e <- matrix(v, ncol = attr(ka, "nrow"))
            attr(e, "type") <- "factor"
            attr(e, "levels") <- tc2
        }
        attr(e, "cellsize") <- attr(ka, "cellsize")
        attr(e, "xll") <- attr(ka, "xll")
        attr(e, "yll") <- attr(ka, "yll")
        class(e) <- "asc"
        ka <- e
    }
    return(ka)
}

