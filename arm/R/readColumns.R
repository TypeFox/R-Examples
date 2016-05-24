read.columns <- function (filename, columns){
    start <- min(columns)
    length <- max(columns) - start + 1
    if (start == 1) {
        return(read.fwf(filename, widths = length))
    }
    else {
        return(read.fwf(filename, widths = c(start - 1, length))[, 2])
    }
}
