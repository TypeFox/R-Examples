where.na <-
function (x) 
{
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    whichna <- seq(along = x)[is.na(x)]
    if (length(whichna) == 0) {
        cat("  No NAs\n")
    }
    else {
        if (is.vector(x)) 
            cat("  Postion(s):", whichna, "\n")
        else {
            rows <- unique(sort(whichna%%dim(x)[[1]]))
            cat("  Row(s):", rows, "\n")
        }
    }
    invisible(whichna)
}
