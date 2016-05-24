


rec <- function(x, slsp=c("remove","missing"))
{
    ## Verifications
    if (!inherits(x, "ltraj"))
        stop("x should be of class \"ltraj\"")

    lif <- infolocs(x)
    if (!is.null(lif)) {
        for (i in 1:length(x)) {
            if (!all(row.names(x[[i]])%in%row.names(lif[[i]]))) {
                x[[i]] <- x[[i]][row.names(x[[i]])%in%row.names(lif[[i]]),,drop=FALSE]
                attr(x[[i]],"infolocs") <- lif[[i]]
            }
            if (!all(row.names(lif[[i]])%in%row.names(x[[i]]))) {
                lif[[i]] <- lif[[i]][row.names(lif[[i]])%in%row.names(x[[i]]),,drop=FALSE]
                attr(x[[i]],"infolocs") <- lif[[i]]
            }
        }
    }

    ## Recomputation
    slsp <- match.arg(slsp)
    if (attr(x,"typeII")) {
        y <- .traj2df(.ltraj2traj(x))
        if (!is.null(lif)) {
            infol <- do.call("rbind", lif)
            al <- as.ltraj(xy=y[,c("x","y")], date=y$date,
                           id=y$id, burst=y$burst, slsp=slsp,
                           typeII=TRUE, infolocs=infol)
        } else {
            al <- as.ltraj(xy=y[,c("x","y")], date=y$date,
                           id=y$id, burst=y$burst, slsp=slsp,
                           typeII=TRUE, infolocs=lif)
        }
        return(al)
    } else {
        attr(x,"typeII") <- TRUE
        y <- .traj2df(.ltraj2traj(x))
        if (!is.null(infolocs(x))) {
            infol <- do.call("rbind", infolocs(x))
            al <- as.ltraj(xy=y[,c("x","y")], id=y$id,
                           burst=y$burst, slsp=slsp,
                           typeII=FALSE, infolocs=infol)
        } else {
            al <- as.ltraj(xy=y[,c("x","y")], id=y$id,
                           burst=y$burst, slsp=slsp,
                           typeII=FALSE)
        }
        return(al)
    }
}
