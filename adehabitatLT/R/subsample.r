### Reechantillonner une trajectoire
subsample <- function(ltraj, dt, nlo=1,
                      units=c("sec", "min", "hour", "day"), ...)
{
    if (!inherits(ltraj,"ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    if ((!is.regular(ltraj))&(attr(ltraj, "typeII")))
        stop("ltraj should be of type I or type II regular")
    if (length(nlo)==1)
        nlo <- rep(nlo, length(ltraj))
    units <- match.arg(units)
    dt <- .convtime(dt, units)

    dtb <- ltraj[[1]]$dt[1]
    if (dt%%dtb!=0)
        stop("dt is not a multiple of the previous time lag")
    la <- dt/dtb
    res <- lapply(1:length(ltraj), function(i) {
        x <- ltraj[[i]]
        infol <- attr(x, "infolocs")
        vec <- rep(1:la, length=nrow(x))
        x <- x[vec==nlo[i],]
        if (!is.null(infol)) {
            infol <- infol[vec==nlo[i],]
            attr(x, "infolocs") <- infol
        }
        return(x)
    })
    class(res) <- c("ltraj","list")
    attr(res,"typeII") <- attr(ltraj,"typeII")
    attr(res,"regular") <- is.regular(res)
    res <- rec(res,...)
    return(res)
}
