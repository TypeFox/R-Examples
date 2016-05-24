"MAX" <- function(
object, lower = 1e-3, upper = 1000, pool = TRUE)
{
    ## Checking class of 'object'
    MAXlist <- object[[11]]$"maxfct"
    if (is.null(MAXlist)) {stop("No method available")}

    ## Retrieving relevant quantities
    indexMat <- object$"indexMat"
    parm <- as.vector(coef(object))
    parmMat <- object$"parmMat"
    strParm <- colnames(parmMat)    
    varMat <- vcov(object, pool = pool)

    ## Calculating ED values    
    ncolIM <- ncol(indexMat)
    indexVec <- 1:ncolIM    
    dimNames <- rep("", ncolIM)
    MAXmat <- matrix(0, ncolIM, 2)

    for (i in indexVec)
    {
        parmInd <- indexMat[,i]
        varCov <- varMat[parmInd, parmInd]
        parmChosen <- parmMat[,i]
        MAXmat[i, ] <- MAXlist(parmChosen, lower, upper)
        dimNames[i] <- strParm[i]
    }

    dimnames(MAXmat) <- list(dimNames, c("Dose", "Response"))
    printCoefmat(MAXmat)
    invisible(MAXmat)    
}
