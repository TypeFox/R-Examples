`calcSensSpec` <-
function(trueMatrix, estMatrix) 
{
    	if(is.matrix(trueMatrix) != TRUE) {
        	stop("Error: ", paste(sQuote("trueMatrix"), sep = ""), 
            	" must be a matrix.")
    	}
    	if (is.matrix(estMatrix) != TRUE) {
        	stop("Error: ", paste(sQuote("estMatrix"), sep = ""), 
            	" must be a matrix.")
    	}
	true <- abs(sign(trueMatrix))
	means <- abs(sign(estMatrix))
	TP <- length(intersect(which(true == 1), which(means == 1)))
	FP <- length(intersect(which(true == 0), which(means == 1)))
	FN <- length(intersect(which(true == 1), which(means == 0)))
	TN <- length(intersect(which(true == 0), which(means == 0)))
	return(list("TP" = TP, "FP" = FP, "FN" = FN, "TN" = TN))
}

