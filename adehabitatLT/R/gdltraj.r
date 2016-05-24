gdltraj <- function(x, min, max,
                    type=c("POSIXct","sec","min","hour",
                    "mday","mon","year","wday","yday"))
{
    ## Verifications
    if (!inherits(x, "ltraj"))
        stop("x should be of class \"ltraj\"")
    if (!attr(x, "typeII"))
        stop("x should be of type II (time recorded)")
    type <- match.arg(type)

    ## gets the traj within the boundaries
    if (type=="POSIXct") {
        x <- lapply(x, function(y) {
            infol <- attr(y, "infolocs")
            if (!is.null(infol))
                infol <- infol[(y$date>min)&(y$date<max),,drop=FALSE]
            y <- y[(y$date>min)&(y$date<max),]
            if (!is.null(infol))
                attr(y, "infolocs") <- infol
            return(y)
        })
    } else {
        x <- lapply(x, function(y) {
            da <- as.POSIXlt(y$date)[[type]]
            infol <- attr(y, "infolocs")
            if (!is.null(infol))
                infol <- infol[(da>=min)&(da<max),,drop=FALSE]
            y <- y[(da>=min)&(da<max),]
            if (!is.null(infol))
                attr(y, "infolocs") <- infol
            return(y)
        })
    }
    if (all(sapply(x,nrow)==0))
        stop("No relocations within the specified interval")
    x[sapply(x, nrow)==0]<-NULL

    ## Output
    class(x) <- c("ltraj", "list")
    attr(x, "typeII") <-  TRUE
    attr(x, "regular") <-  is.regular(x)
    x <- rec(x)
    return(x)
  }
