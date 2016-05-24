epi.asc <- function(dat, file, xllcorner, yllcorner, cellsize, na = -9999) 
{
    id <- is.na(dat)
    dat[id] <- na
    ncols <- dim(dat)[2]
    nrows <- dim(dat)[1]
    h.ncol <- paste("ncols", nrows)
    h.nrow <- paste("nrows", ncols)
    # h.ncol <- paste("ncols", ncols)
    # h.nrow <- paste("nrows", nrows)
    h.xllcorner <- paste("xllcorner", xllcorner)
    h.yllcorner <- paste("yllcorner", yllcorner)
    h.cellsize <- paste("cellsize", cellsize)
    h.nodata <- paste("nodata_value", na)
    header <- rbind(h.ncol, h.nrow, h.xllcorner, h.yllcorner, h.cellsize, h.nodata)
    write.table(header, file = file, append = FALSE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

    rval <- as.matrix(dat)
    rval <- matrix(rval, nrow = ncols, byrow = TRUE)
    rval <- rval[ncols:1,]

    write.table(rval, file = file, append = TRUE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
}
