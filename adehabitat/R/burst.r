burst <- function(ltraj)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    return(unlist(lapply(ltraj, function(x) attr(x,"burst"))))
}##OK

"burst<-" <- function(ltraj, value)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    res <- lapply(1:length(ltraj), function(i) {
        x <- ltraj[[i]]
        attr(x,"burst") <- value[i]
        return(x)
    })
    class(res) <- c("ltraj","list")
    attr(res, "typeII") <- attr(ltraj, "typeII")
    attr(res, "regular") <- is.regular(res)
    return(res)
}##OK


id <- function(ltraj)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    return(unlist(lapply(ltraj, function(x) attr(x,"id"))))
}##OK

"id<-" <- function(ltraj, value)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    res <- lapply(1:length(ltraj), function(i) {
        x <- ltraj[[i]]
        attr(x,"id") <- value[i]
        return(x)
    })
    class(res) <- c("ltraj","list")
    attr(res, "typeII") <- attr(ltraj, "typeII")
    attr(res, "regular") <- is.regular(res)
    return(res)
}##OK
