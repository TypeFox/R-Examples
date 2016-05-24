#' @name dunnetts.format
#' @title Format Dunnett's Test Output
#' @param summary  Summary object from dunnett's test for formatting
#' @param test   Name of the test to be formatted
#' @param levels   Dose factor levels defining rows in formatted output
#' @keywords internal

dunnetts.format <- function (summary, test, levels) {

    coeffchar <- as.character (summary$test$coefficients)
    sigmachar <- as.character (summary$test$sigma)
    tchar <- as.character(summary$test$tstat)
    pchar <- summary$test$pvalues
   
    
    outputlist <- list()
    
    # Format output into matrices, 1 for the title, 1 for the results
    # Title
    title <- matrix(nrow=1, ncol=1)
    title[1,1] <- test
    
    # Results
    output_rows <- length(levels) + 1
    
    outputmatrix <- matrix(nrow=output_rows, ncol=5)
	outputmatrix[1,1] <- "Dose vs Zero "
    outputmatrix[1,2] <- "Estimate "
    outputmatrix[1,3] <- "Std Error "
    outputmatrix[1,4] <- "t value "
    outputmatrix[1,5] <- "p value"

	iterations <- length(levels)
	for (i in 1:iterations)
	{
    outputmatrix[i+1,1] <- strtrim(levels[i+1], 6)
    outputmatrix[i+1,2] <- strtrim(coeffchar[i], 6)
    outputmatrix[i+1,3] <- strtrim(sigmachar[i], 6)
	outputmatrix[i+1,4] <- strtrim(tchar[i], 6)
 	outputmatrix[i+1,5] <- round(pchar[i], 10)
	}
    outputlist[[1]] <- title
    outputlist[[2]] <- outputmatrix
    return(outputlist)
}