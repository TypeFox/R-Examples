commatFct <- function(object, compMatch)
{
    parmMat <- object$parmMat

    if (!is.null(compMatch))
    {
        return(parmMat[, (colnames(parmMat) %in% c(compMatch[1], compMatch[2])), drop = FALSE ])
    } else {
        parmMat
    }
}