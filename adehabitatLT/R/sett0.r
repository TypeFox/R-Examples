## Pour poser le t0 en fonction dune date de reference
sett0 <- function(ltraj, date.ref, dt,
                  correction.xy=c("none", "cs"), tol=dt/10,
                  units=c("sec", "min", "hour", "day"), ...)
{
    if (!inherits(ltraj,"ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    if (!attr(ltraj, "typeII"))
        stop("ltraj should be of type II (time recorded)")
    if (inherits(date.ref,"POSIXlt"))
        date.ref <- as.POSIXct(date.ref)
    if (is.integer(date.ref)) {
        class(date.ref) <- c("POSIXct","POSIXt")
        attr(date.ref, "tzone") <- attr(ltraj[[1]]$date, "tzone")
    }
    if (!inherits(date.ref,"POSIXct"))
        stop("date.ref should be of class \"POSIXct\"")

    units <- match.arg(units)
    dt <- .convtime(dt, units)
    tol <- .convtime(tol, units)
    if (length(date.ref)==1)
        date.ref <- rep(date.ref, length(ltraj))

    correction.xy <- match.arg(correction.xy)
    res <- lapply(1:length(ltraj), function(oo) {
        x <- ltraj[[oo]]
        infol <- attr(x, "infolocs")
        date.reft <- date.ref[oo]
        dc <- x$date
        da <- as.numeric(x$da) - as.numeric(date.reft)
        x$date <- round(da/dt,0)*dt + as.numeric(date.reft)
        if (any(abs(as.numeric(dc) - as.numeric(x$date)) > tol))
            stop("ltraj contains irregular data (time lag > or < tol)")
        class(x$date) <- c("POSIXct","POSIXt")
        attr(x$date, "tzone") <- attr(ltraj[[1]]$date, "tzone")
        if (correction.xy=="cs") {
            rr <- .corrXY(x$x, x$y,
                          as.numeric(dc),
                          as.numeric(x$date))
            x$x <- rr$x
            x$y <- rr$y
        }
        attr(x, "infolocs") <- infol
        return(x)
    })
    class(res) <- c("ltraj", "list")
    attr(res, "typeII") <- TRUE
    attr(res, "regular") <- TRUE
    res <- rec(res,...)
    return(res)
}
