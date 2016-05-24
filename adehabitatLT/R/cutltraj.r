
cutltraj <- function(ltraj, criterion, value.NA = FALSE,
                     nextr = TRUE, ...)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    infol <- infolocs(ltraj)
    res <- list()
    k <- 1
    for (i in 1:length(ltraj)) {
        if (!is.null(infol)) {
            att <- attributes(ltraj[[i]])
            x <- cbind(ltraj[[i]], infol[[i]])
            attr(x, "id") <- att$id
            attr(x, "burst") <- att$burst
        } else {
            x <- ltraj[[i]]
        }
        ex <- parse(text = criterion)
        coin <- eval(ex, envir = x)
        coin[is.na(coin)] <- value.NA
        if (nextr) {
            kwa <- c(1, cumsum(as.numeric(coin))+1)
            kwa <- kwa[1:(length(kwa)-1)]
        } else {
            kwa <- cumsum(as.numeric(coin))
            x <- x[!coin,]
            kwa <- kwa[!coin]
        }
        mkk <- nchar(max(kwa))
        kwa <- sapply(kwa, function(hh) {
            nc <- nchar(hh)
            if (mkk-nc>0) {
                return(paste(c(rep("0",mkk-nc),hh), collapse=""))
            } else {
                return(as.character(hh))
            }})
        bu <- factor(paste(attr(x,"burst"),kwa,sep="."))
        id <- factor(rep(attr(x,"id"), nrow(x)))

        if (is.null(infol)) {
            res[[i]] <- as.ltraj(id=id, xy=x[,c("x","y")],
                                 date=x$date, burst=bu,
                                 typeII=attr(ltraj,"typeII"),
                                 infolocs=infol)
        } else {
            inf <- x[,(names(x)%in%names(infol[[i]])), drop=FALSE]
            res[[i]] <- as.ltraj(id=id, xy=x[,c("x","y")],
                                 date=x$date, burst=bu,
                                 typeII=attr(ltraj,"typeII"),
                                 infolocs=inf)
        }
    }
    if (length(res)>1) {
        res <- do.call("c.ltraj", res)
    } else {
        res <- res[[1]]
    }
    rrr <- unlist(lapply(res,nrow))
    resb <- res[rrr>=3]
    if (length(res)!=length(resb))
        warning(paste("At least 3 relocations are needed for a burst\n",
                      sum(rrr[rrr<3]), "relocations have been deleted"))
    resb <- rec(resb,...)
    return(resb)
}

bindltraj <- function(ltraj, ...)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    traj <- .ltraj2traj(ltraj)
    traj$burst <- traj$id
    return(.traj2ltraj(traj, ...))
}
