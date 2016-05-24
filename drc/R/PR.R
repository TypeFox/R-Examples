"PR" <- function(object, xVec, ...)
{
    lenXV <- length(xVec)
    
    curveId <- as.character(unique(object$data[, 3]))
    lenCI <- length(curveId)
    
    if (lenCI > 1)
    {
        retMat <- predict(object, data.frame(xVec, rep(curveId, rep(lenXV, lenCI))), se.fit = TRUE, ...)
        rownames(retMat) <- paste(rep(curveId, rep(lenXV, lenCI)), rep(as.character(xVec), lenCI), sep = ":")
    } else {
        retMat <- predict(object, data.frame(xVec))
        if (is.matrix(retMat))
        {
            rownames(retMat) <- rep(as.character(xVec), lenCI)
        } else {
            names(retMat) <- rep(as.character(xVec), lenCI)
        }
    }
    
    return(retMat)
}
