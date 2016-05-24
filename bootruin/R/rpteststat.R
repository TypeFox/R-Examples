rpteststat <- function(x, x.null, se){
    stopifnot(is.numeric(x), is.numeric(x.null), is.numeric(se), is.vector(x.null))

    if (is.vector(x) && is.vector(se)) {
        stopifnot((nc <- length(x)) == length(se))
        nr <- length(x.null)
        Recall(x      = matrix(x,  ncol = nc, nrow = nr, byrow = TRUE),
               x.null = x.null,
               se     = matrix(se, ncol = nc, nrow = nr, byrow = TRUE))
    } else {
        stopifnot(is.matrix(x), is.matrix(se), dim(x) == dim(se))
        sweep(x, 1L, x.null) / se
    }
}
