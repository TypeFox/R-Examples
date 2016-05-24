FBN.valueCenter <-
function(inputData = NULL, normalizingValue = NULL, nominalValueCN = 2, logScale = FALSE){	
	if(is.null(inputData)){
        cat("WARNING: FBN.valueCenter -> Please input a valid inputData\n") 
        return(NULL)
    }
    if(is.null(normalizingValue) ){
        cat("WARNING: FBN.valueCenter -> Please input a valid normalizingValue\n") 
        return(NULL)
    }
    if(logScale)
    	normalizedData = inputData + log2(nominalValueCN) - log2(normalizingValue) 
    else
		normalizedData = inputData * nominalValueCN / normalizingValue
	return(normalizedData)
}

