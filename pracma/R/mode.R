##
##  m o d e . R
##


Mode <- function(x) {
    if (is.matrix(x))
        x <- c(x)

    if (is.numeric(x)) {
        x   <- sort(x)
        tbl <- table(x)
        n   <- which.max(tbl)
        xm  <- as.numeric(names(tbl)[n])
    } else if (is.complex(x)) {
        x   <- x[order(abs(x))]
        tbl <- table(x)
        n   <- which.max(tbl)
        xm  <- as.complex(names(tbl)[n])
    } else if (is.factor(x)) {
        tbl <- table(x)
        n   <- which.max(tbl)
        xm  <- names(tbl)[n]        
    } else
        xm <- NA

    return(xm)
}
