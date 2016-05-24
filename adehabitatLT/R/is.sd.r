is.sd <- function(ltraj)
{
    if (!inherits(ltraj,"ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    if (!is.regular(ltraj))
        stop("ltraj should be regular")
    len <- unlist(lapply(ltraj, nrow))
    len <- all(len==len[1])
    return(len)
}

sd2df <- function(ltraj, what)
{
    if (!inherits(ltraj,"ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    if (!is.sd(ltraj))
        stop("ltraj should contain regular bursts (same time lag), with the same duration")

    res <- lapply(ltraj, function(x) {
        if (!is.null(attr(x,"infolocs")))
            x <- cbind(x, attr(x,"infolocs"))
        x[what]
    })
    res <- do.call("data.frame",res)
    names(res) <- burst(ltraj)
    return(res)
}
