read.ifile <- function(filename, skip = 1, col.names,
    sep = ",", ...) {

    out <- read.table(filename, skip = skip, sep = sep, ...)

    if(!missing(col.names))
        colnames(out) <- col.names
    else {
        colnames(out) <- gsub("X.", "", colnames(out))
        colnames(out) <- tolower(colnames(out))
    }

    out <- as.ifile(out)

    return(out)
}
