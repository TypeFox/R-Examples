`trans.comprisk` <- function(K, names)
{
    tmat <- matrix(NA, K+1, K+1)
    tmat[1, 2:(K+1)] <- 1:K
    if (missing(names))
        names <- c("eventfree",paste("cause",1:K,sep=""))
    else {
        if (length(names)==K) names <- c("eventfree",names)
        else {
            if (length(names) != K+1)
                stop("incorrect length of \"names\" argument")
        }
    }
    dimnames(tmat) <- list(from=names, to=names)
    return(tmat)
}

`trans.illdeath` <- function(names)
{
    tmat <- matrix(NA, 3, 3)
    tmat[1, 2:3] <- 1:2
    tmat[2, 3] <- 3
    if (missing(names))
        names <- c("healthy","illness","death")
    else {
        if (length(names)!=3) stop("incorrect length of \"names\" argument")
    }
    dimnames(tmat) <- list(from=names, to=names)
    return(tmat)
}
