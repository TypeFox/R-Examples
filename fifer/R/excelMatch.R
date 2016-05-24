##' Generate a sequence of column labels to match Excel's naming
##'
##' Excel columns are labeled with letters (e.g., A, B, C, ... AA, AB, AC, etc). Given an integer
##' (n), this function will return labels starting from A until the nth column. See examples.
##' @title Generate Excel column labels
##' @param n an integer that indicates how many named columns the user wishes to obtain
##' @return a vector of strings of length n
##' @author Dustin Fife
##' @seealso \code{\link{excelMatch}}
##' @examples
##' excelCols(11)
##' @export
excelCols=function(n){
	vec1=LETTERS
	if(n<=26){
		res=vec1[seq_len(n)]
	} else if(n>26&n<=702){
		res=c(vec1,apply(expand.grid(vec1,vec1)[,2:1],1,paste,collapse=""))[1:n]
	} else if(n>702&n<=18278){
		res=c(vec1,apply(expand.grid(vec1,vec1)[,2:1],1,paste,collapse=""),apply(expand.grid(vec1,vec1,vec1)[,3:1],1,paste,collapse=""))[1:n]
	} else{
		res=NA###forgot
	}
	res
}

##' Obtain column number or variable name from Excel named Columns
##'
##' Excel columns are labeled with letters (e.g., A, B, C, ... AA, AB, AC, etc). This function
##' accepts a string (e.g., "AAC") and returns either a number that indicates where that string falls
##' in the sequence of excel named columns, or it returns the variable name corresponding to that column number.
##' 
##' @title Obtain column number or variable name from Excel named Columns
##' @param ... An Excel-like string consisting of all capital letters (e.g., "AAQ", "BZ", "RQ", "S", etc.)
##' @param n the number of columns of the data of interest
##' @param names the column names of the data of interest
##' @return either a variable name or column number
##' @seealso \code{\link{excelCols}}
##' @examples
##' fake.names = paste("Variable", 1:1000, sep="")
##' # find the Variable name corresponding to AC
##' excelMatch("AC", names=fake.names)
##' # find the column number instead
##' excelMatch("AC", n=1000)
##' @author Dustin Fife
##' @export
excelMatch = function(..., n=NULL, names=NULL){

	#### check for errors
    if (is.null(n) & is.null(names)){
    	stop("You must supply either an n value (integer) or a vector of names(names)")	
    } else if (!is.null(names)){
	    n = length(names)    	
    }    

	#### generate excel columns
	excelColumns = excelCols(n)
	
	#### do the search
	findit = which(excelColumns %in% unlist(list(...)))

	if (!is.null(names)){
		return(names[findit])
	} else {
		return(findit)
	}
}