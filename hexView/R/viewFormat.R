
# User functions for creating and viewing rawFormat's

readFormat <- function(file,
                       format,
                       width=NULL, offset=0, 
                       machine="hex",
                       flatten=TRUE) {
    if (!is.memFormat(format))
        stop("Invalid format")

    fileSize <- fileSize(file)

    infile <- file(file, "rb")
    on.exit(close(infile))

    if (offset > 0)
        seek(infile, offset)

    blocks <- lapply(format, readBlock, infile)
    fileFormat <- list(blocks=blocks)
    if (flatten) {
        fileFormat$blocks <- flattenFormat(fileFormat$blocks)
        class(fileFormat) <- c("flatRawFormat", "rawFormat")
    } else {
        class(fileFormat) <- "rawFormat"
    }
    fileFormat$offset <- offset
    fileFormat$nbytes <- fileSize - offset
    fileFormat
}


viewFormat <- function(..., page=FALSE) {
    print(readFormat(...), page=page)
}
