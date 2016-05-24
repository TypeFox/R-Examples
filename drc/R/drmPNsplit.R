"drmPNsplit" <- 
function(parmVec, sep) 
{
    lenPV <- length(parmVec)
    parmVecA <- rep(0, lenPV)
    parmVecB <- rep(0, lenPV)

    splitList <- strsplit(parmVec, sep, fixed = TRUE)
    for (i in 1:lenPV)
    {
        parmVecA[i] <- splitList[[i]][1]
        
        lenSL <- length(splitList[[i]])
        parmVecB[i] <- paste(splitList[[i]][2:lenSL], collapse = "")  # 'paste' is needed in case several ":" occur
    }
    return(list(parmVec, parmVecA, parmVecB))
}
