permuteGenotype <-
function (dataLine, fromCol) 
{
    # if (fromCol<0 || fromCol != round(fromCol)) {
        # stop("'fromCol' must be a positive integer.")
    # }
    colNum = ncol(dataLine)
    randIndex = sample(colNum - fromCol+1)
	newDataLine=cbind(dataLine[,1:(fromCol-1),drop=FALSE],dataLine[,randIndex+fromCol-1,drop=FALSE])
	return(newDataLine)
}
