EtaSq<-function (x) 
{
    anovaResults <- summary.aov(x)[[1]]
    anovaResultsNames <- rownames(anovaResults)
    SS <- anovaResults[,2] #SS effects and residuals
    k <- length(SS) - 1  # Number of factors 
    ssResid <- SS[k + 1]  # Sum of Squares Residual
    ssTot <- sum(SS)  # Sum of Squares Total
    SS <- SS[1:k] # takes only the effect SS
    anovaResultsNames <- anovaResultsNames[1:k]
    etaSquared <- SS/ssTot # Should be the same as R^2 values
    partialEtaSquared <- SS/(SS + ssResid)
    res <- cbind(etaSquared, partialEtaSquared)
    colnames(res) <- c("Eta^2", "Partial Eta^2")
    rownames(res) <- anovaResultsNames
    return(res)
}
