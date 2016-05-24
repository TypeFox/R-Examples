"Rsq" <- function(object)
{
    response <- object$data[,2]
    curve <- object$data[,4]
    uniCurve <- unique(curve)
    lenUC <- length(uniCurve)

    numerator <- tapply( residuals(object)^2, curve, sum)  # residual ss
    denominator <- tapply( (response - mean(response))^2, curve, sum)  # total SS

    totnum <- sum(residuals(object)^2)
    totden <- sum((response - mean(response))^2)

    if (lenUC==1)
    {
        hText <- "\nR-square value\n"
        rsq <- matrix(c(1 - numerator/denominator), 1, 1)
        rownames(rsq) <- "" 
    } else {
        hText <- "\nR-square values\n"
        rsq <- matrix(c(1 - numerator/denominator, 1-totnum/totden), lenUC+1, 1)
        rownames(rsq) <- c(as.character(uniCurve), "Total") 
    }
    colnames(rsq) <- ""
    
    cat(hText)
    printCoefmat(rsq)
    invisible(rsq)
}
