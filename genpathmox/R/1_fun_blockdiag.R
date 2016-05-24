#' @title Building the Bart matrix
#' @details
#' Internal function. \code{blockdiag} is called by \code{build.block}.
#' @param x list of matrix used to build the Bart matrix
#' @param \dots Further arguments passed on to \code{\link{blockdiag}}.
#' @return Bart matrix
#' @import diagram methods mice plspm quantreg
#' @keywords internal
#' @export

blockdiag <- function (x, ...)
{
    if (!is.list(x)) 
        x <- list(x)
    args <- list(...)
    if (length(args) > 0) 
        args <- c(x, args)
    else args <- x
    idx <- which(!sapply(args, is.matrix))
    if (length(idx) > 0) 
        for (i in idx) args[[i]] <- as.matrix(args[[i]])
    if (length(args) == 1) 
        return(args[[1]])
    nr <- sapply(args, nrow)
    nc <- sapply(args, ncol)
    cumnc <- cumsum(nc)
    NR <- sum(nr)
    NC <- sum(nc)
    rowfun <- function(m, zbefore, zafter) {
        cbind(matrix(0, ncol = zbefore, nrow = nrow(m)), m, matrix(0, 
            ncol = zafter, nrow = nrow(m)))
    }
    ret <- rowfun(args[[1]], 0, NC - ncol(args[[1]]))
    for (i in 2:length(args)) {
        ret <- rbind(ret, rowfun(args[[i]], cumnc[i - 1], NC - 
            cumnc[i]))
    }
    ret
}
