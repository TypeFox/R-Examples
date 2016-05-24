matchVertCons <-
function (fullhull)
{
    if (!is.null(fullhull[["affectedConstraints"]]))
        return(fullhull)
    gmat <- t(fullhull$a %*% t(fullhull$vertices) - fullhull$b)
    testmat <- (gmat == 0)
    dim(testmat) <- dim(gmat)
    affectedConstraints <- apply(testmat, 1, which)
    return(c(fullhull, list(affectedConstraints = affectedConstraints)))
}
