"drmParNames" <- 
function(numNames, parNames, collapseList2, repStr1 = "factor(pmodels[, i])", repStr2 = "factor(assayNo)")
{
    ## Retrieving names for parameters
    parmVecList <- list()
    for (i in 1:numNames)
    {
        colNames1 <- colnames(collapseList2[[i]])
        if (is.null(colNames1)) 
        {
            parmVecList[[i]] <- paste(parNames[i], "(Intercept)", sep = ":")
        } else { 
            parmVecList[[i]] <- paste(parNames[i], colNames1, sep = ":")
        }        
        parmVecList[[i]] <- (parmVecList[[i]])[1:ncol(collapseList2[[i]])]  # min(maxParm[i], length(colNames1))]
    }
    parmVec <- unlist(parmVecList)
        
    parmVec2 <- parmVec
#    print(parmVec2)
    for (i in 1:length(parmVec))
    {
        pos <- regexpr(repStr1, parmVec[i], fixed = TRUE)
        if (pos > 0) 
        {
            parmVec2[i] <- paste(substring(parmVec[i], 1, pos-1), substring(parmVec[i], pos + 20), sep = "")
        }
            
        pos <- regexpr(repStr2, parmVec[i], fixed = TRUE)
        if (pos > 0) 
        {
            parmVec2[i] <- paste(substring(parmVec[i], 1, pos-1), substring(parmVec[i], pos + 15), sep = "")
        }
    }
    
#    print(parmVec2)
    return(drmPNsplit(parmVec2, ":"))
}
