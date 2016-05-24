ld <- function(ltraj)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class ltraj")
    iidd <- rep(id(ltraj), sapply(ltraj, nrow))
    bur <- rep(burst(ltraj), sapply(ltraj, nrow))
    inf <- infolocs(ltraj)
    tr <- do.call("rbind", ltraj)
    if (!is.null(inf))
        return(data.frame(tr, id=iidd, burst=bur, do.call("rbind",inf)))
    return(data.frame(tr, id=iidd, burst=bur))
}

dl <- function(x)
{
    if (!inherits(x, "data.frame"))
        stop("x should be of class data.frame")
    nareq <- c("x","y","date")
    sapply(nareq, function(j) {
        if (!any(names(x)==j))
            stop(paste("No variable named", j))
    })
    if (!any(names(x)=="id")) {
        idd <- bur <- rep(paste(sample(letters, 5), collapse=""), nrow(x))
    } else {
        idd <- x$id
    }
    if (!any(names(x)=="burst")) {
        bur <- idd
    } else {
        bur <- x$burst
    }

    xy <- x[,c("x","y")]
    dat <- x$date
    type2 <- TRUE
    if (!inherits(dat, "POSIXct"))
        type2 <- FALSE
    pasinf <- c("x","y","date","dx", "dy", "dist", "dt", "R2n", "abs.angle", "rel.angle", "id","burst")
    x <- x[,!(names(x)%in%pasinf), drop=FALSE]
    if (ncol(x)>0) {
        inf <- x
    } else {
        inf <- NULL
    }
    as.ltraj(xy, dat, idd, bur, type2, infolocs=inf)
}
