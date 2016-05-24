fcrosRead <-
function(filename) {
    xdata <- read.table(filename, sep="\t", header = TRUE)
    return(xdata)
}
