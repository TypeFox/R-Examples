rpdataconv <- function(x) {
    stopifnot(all(vapply(x, is.numeric, logical(1L))))
    if (is.vector(x) && !is.list(x)) {
        return(matrix(x, ncol = 1L))
    } else {
        if (is.matrix(x)) {
            return(x)
        } else {
            x       <- lapply(x, as.vector)
            x.len   <- vapply(x, length, integer(1L))
            max.len <- max(x.len)
            if(max.len == 1L){
                warning("All entries of x have length 1.")
                return(matrix(unlist(unname(x)), nrow = 1L))
            } else {
                x.fill <- lapply(max.len - x.len, rep.int, x = NA_real_)
                x      <- mapply(FUN = c, x, x.fill)
                return(unname(as.matrix(as.data.frame(x))))
            }
        }
    }
}
