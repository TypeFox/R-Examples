## Copied from mixstock
blockdiag <- function (...) 
{
    args <- list(...)
    nc <- sapply(args, ncol)
    cumnc <- cumsum(nc)
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
