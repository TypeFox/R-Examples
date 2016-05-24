"print.kasc" <- function(x, ...)
{
    ## Verifications
    if (!inherits(x, "kasc")) stop("Non convenient data")

    ## The output
    cat("Raster map of class \"kasc\":\n")
    cat("Cell size: ", attr(x, "cellsize"), "\n")
    cat("Number of rows: ", attr(x, "nrow"), "\n")
    cat("Number of columns: ", attr(x, "ncol"), "\n\n")

    cat("Variables measured:\n")
    n<-names(x)
    for (i in 1:length(n)) {
        if (is.factor(x[[i]])) {
            typ<-"factor"
        } else {
            typ<-"numeric"
        }
        cat(paste(i, ". ", n[i], ": ", typ, "\n", sep=""))
    }
    cat("\n")
}

