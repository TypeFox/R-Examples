##' Given a dataset, \code{missing.vals} will report on the number
##' of missing observations per variable
##'
##' @title Summary of Missing Data
##' @param dataset a matrix or data frame (persons in rows, variables in columns).
##' @return a matrix with information about what is missing.
##' @export
##' @aliases missingVals missing.Vals missingvals
##' @author Dustin Fife
missing.vals = function(dataset){
	ms = sapply(dataset, function(x) sum(is.na(x)))
	ms = ms[ms>0]
	ord = order(ms, decreasing=TRUE)
	m = data.frame(matrix(ms[ord], nrow=length(ms), ncol=1))
	names(m) = "Number Missing"
	row.names(m) = names(ms)[ord]
	if (nrow(m)<1){
		cat("There are no missing values in this dataset.")
	} else {
		mx = max(m[,1])
		mn = min(m[,1])
		return(m)	
	}
}	
