`fill.na` <-
function (x)
{
    if (NCOL(x) == 1) {
        x <- as.matrix(x)
        onedim <- TRUE
    } else onedim <- FALSE
    if (any(is.na(x[1,])))
        stop("cannot replace 'NA's in first place")
    out <- x
    for (cols in 1:ncol(x)) {
        for (rows in 1:nrow(x)) {
            if (is.na(x[rows, cols])) {
                rows2 <- rows - 1
                out[rows, cols] <- out[rows2, cols]
            }
        }
    }
    if (onedim)
        out <- array(out)
    return(out)
}

