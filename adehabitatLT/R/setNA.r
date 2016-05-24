setNA <- function(ltraj, date.ref, dt, tol=dt/10,
                  units=c("sec", "min", "hour", "day"), ...)
{
    if (!inherits(ltraj,"ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    if (!attr(ltraj, "typeII"))
        stop("ltraj should be of type II (time recorded)")
    if (is.numeric(date.ref)) {
        class(date.ref) <- c("POSIXct","POSIXt")
        attr(date.ref, "tzone") <- attr(ltraj[[1]]$date, "tzone")
    }
    if (inherits(date.ref,"POSIXlt"))
        date.ref <- as.POSIXct(date.ref)
    if (!inherits(date.ref,"POSIXct"))
        stop("date.ref should be of class \"POSIXct\"")

    units <- match.arg(units)
    dt <- .convtime(dt, units)
    tol <- .convtime(tol, units)
    if (length(date.ref)==1)
        date.ref <- rep(date.ref, length(ltraj))

    res <- lapply(1:length(ltraj), function(oo) {
         x <- ltraj[[oo]]
         infol <- attr(x, "infolocs")
         date.refp <- date.ref[oo]
         dc <- x$date
         da <- as.numeric(x$da) - as.numeric(date.refp)
         glou <- round(da/dt,0)*dt + as.numeric(date.refp)
         if (any(abs(as.numeric(dc) - as.numeric(glou)) > tol))
            stop("ltraj contains irregular data (time lag > or < tol)")
         laou <- as.integer(round(da/dt,0))
         mlaou <- min(laou)
         laou <- laou-mlaou+1
         xx <- rep(NA, max(laou))
         yy <- rep(NA, max(laou))
         da <- (((1:max(laou))-1)+mlaou)*dt + as.numeric(date.refp)
         if (!is.null(infol)) {
             infol <- do.call("data.frame", lapply(infol, function(y) {
                 ll <- rep(NA, max(laou))
                 ll[laou] <- y
                 return(ll)
             }))
         }
         xx[laou] <- x$x
         yy[laou] <- x$y
         da[laou] <- x$date
         class(da) <- c("POSIXct","POSIXt")
         attr(da, "tzone") <- attr(ltraj[[1]]$date, "tzone")
         return(as.ltraj(data.frame(xx,yy), da, id=attr(x,"id"),
                         burst = attr(x,"burst"), typeII=TRUE,
                         infolocs=infol, ...))
     })
    return(do.call("c.ltraj",res))
}
