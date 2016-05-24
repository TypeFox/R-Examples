#' @name dunnettst3
#' @title Dunnett's T3 Test
#' @param targetcolumn  Character string, name of response column to be tested
#' @param alternative Character string(s) specifying the direction of the alternative hypothesis.
#' @param alpha  Significance level (numeric) to be used.
#' @param control  Not relevant for this function
#' @param tot.obs  Not relevant for this function
#' @param label  Label of the alternative direction for output.
#' @param data  Input dataframe.
#' @keywords internal

dunnettst3 <- function (targetcolumn, alternative, alpha, control, tot.obs, label, data) {
  # Determine alpha level from label passed in from noel function:
  if (length(grep("One-Tailed", label)) > 0) {
    alpha <- alpha*2
    label <- "One-Tailed"}

	dunnett_t3 <- DTK::DTK.test(x=data[,targetcolumn], f=data$dose_fac, a=alpha)
    levels <- levels(data$dose_fac)
    redlevels <- levels[levels != "0"]
    fullresultmat <- dunnett_t3[[2]]
    redresultmat <- fullresultmat[1:length(redlevels), 1:3] 
    
    outputmatrix <- matrix(nrow=(length(redlevels)+3), ncol=4)
    outputmatrix[1,1] <- paste(label, "Dunnett's T3 Test", sep = " ")
    outputmatrix[3,1] <- "Dose Levels"
    outputmatrix[3,2] <- "Difference"
    outputmatrix[3,3] <- "Lower CI"
    outputmatrix[3,4] <- "Upper CI"
    
    for (i in 1:length(redlevels)) {
        
        outputmatrix[i+3,1] <- strtrim(levels[[i+1]], 6)
        outputmatrix[i+3,2] <- round(redresultmat[i,1], digits = 4)
        outputmatrix[i+3,3] <- round(redresultmat[i,2], digits = 4)
        outputmatrix[i+3,4] <- round(redresultmat[i,3], digits = 4)
    }
    
    outputlist <- list()
    outputlist[[1]] <- outputmatrix
	return(outputlist)
}
