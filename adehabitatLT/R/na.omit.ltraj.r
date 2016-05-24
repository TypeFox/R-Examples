
na.omit.ltraj <- function(object, ...)
{
    if (!inherits(object, "ltraj"))
        stop("ltraj should be of class ltraj")

    info <- infolocs(object)
    for (i in 1:length(object)) {
        x <- object[[i]]
        if (!is.null(info))
            info[[i]] <- info[[i]][(!is.na(x[,1]))&(!is.na(x[,2])),,drop=FALSE]
        object[[i]] <- x[(!is.na(x[,1]))&(!is.na(x[,2])),]
    }
    if (!is.null(info))
        infolocs(object) <- info
    return(rec(object))
}
