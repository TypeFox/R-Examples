
cutltraj <- function(ltraj, criterion, value.NA = FALSE,
                     nextr = TRUE, ...)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    res <- list()
    k <- 1
    for (i in 1:length(ltraj)) {
        x <- ltraj[[i]]
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
        res[[i]] <- as.ltraj(id=id, xy=x[,c("x","y")],
                             date=x$date, burst=bu,
                             typeII=attr(ltraj,"typeII"))
    }
    res <- do.call("c.ltraj", res)
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
    traj <- ltraj2traj(ltraj)
    traj$burst <- traj$id
    return(traj2ltraj(traj, ...))
}
