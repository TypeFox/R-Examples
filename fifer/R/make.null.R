##' Given a starting string and an ending string, \code{r} will extract all columns between the first 
##' and the last string and return either a vector of strings or a vector of integers.
##' See examples. 
##'
##' @title Extract variable column indices
##' @param stringA the name of the first variable you wish to extract the column for in sequence
##' @param stringB the name of the last variable you wish to extract the column for in sequence
##' @param data.names the names of the dataset
##' @param names Should the names of the variables be returned? Or the column indices? Defaults to FALSE 
##' (meaning the column indices will be returned).
##' @return a vector of numbers corresponding to the column names
##' @seealso \code{\link{make.null}}, \code{\link{get.cols}}
##' @author Dustin Fife
##' @examples
##' var.names = LETTERS[1:20]
##' r("C", "F", var.names)
##' @export
r = function(stringA, stringB, data.names, names=F){
	
	#### make sure those names exist
	start = which(data.names==stringA)
	end = which(data.names==stringB)
	
	##### give warning if they can't find these names
	if (length(start)<1){
		warning(paste("\nCould not find the variable", stringA, "in the names"))
		stop()
	} 
	if (length(end)<1){
		warning(paste("\nCould not find the variable", stringB, "in the names"))
		stop()
	} 
	if (names){
		return(data.names[start:end])
	} else {
		return(start:end)
	}
}


##' Given a vector of strings and another reference vector, \code{get.cols} will search the reference
##' vector for matches and return the column indices where each string was located. See examples. 
##'
##' @title Extract column index
##' @param string the name of a variable (string)
##' @param names the names of the variables for which you wish to extract the column
##' @return an integer corresponding to the column index(indices)
##' @author Dustin Fife
##' @export 
##' @seealso \code{\link{make.null}}, \code{\link{r}}
##' @examples
##' names = LETTERS[1:10]
##' get.cols(c("A", "B", "D"), LETTERS[1:10])  ### works
##' #get.cols(c("A", "B", "Z"), LETTERS[1:10])  ### doesn't work
get.cols = function(string, names){
	
	#### first check to make sure they're all there
	check = which(!(string %in% names))
	if (length(check)>0){
		msg = "The following variables could not be located in the names string:\n"
		msg2 = paste(string[check], collapse=",")
		warning(paste(msg, msg2))
		stop()
		
	}
	a = which(names %in% string)
	return(a)
}



##' Given a set of variable names (or integers), remove (or keep) these columns from a dataset
##'
##' @title Drop or keep variables in a dataset
##' @export
##' @param ... objects (either vectors of strings, vectors of numbers, or both)
##' @param data the dataset you're trying to eliminate (or keep) variables from
##' @param keep should the variables be kept (keep=T) or dropped (keep=F)?
##' @seealso \code{\link{subset}}. 
##' @return a dataset, containing the variables of interest only
##' @author Dustin Fife
##' @seealso \code{\link{get.cols}}, \code{\link{r}}
##' @examples
##' data = data.frame(matrix(rnorm(100), ncol=5))
##' names(data) = LETTERS[1:5]
##' 
##' #### extract only the classification
##' data(iris)
##' new.data = make.null(c("Sepal.Length"), 
#'        r("Sepal.Width", "Petal.Width", names(iris)), data=iris)
##'
##' #### extract all but the classification
##' new.data2 = make.null(c("Sepal.Length"), 
##'    r("Sepal.Width", "Petal.Width", names(iris)), data=iris, keep=TRUE)
make.null = function(..., data, keep=FALSE){
	
	#### get arguments
	args = list(...)
	
	#### extract names
	nms = names(data)
	
	#### figure out which is numeric/character
	modes = sapply(list(...), mode)
	char.ones = which(modes=="character")
	numb.ones = which(modes=="numeric")	

	#### convert names to column numbers
	colnumbs = sapply(args[char.ones], get.cols, names=nms)
	
	#### put numbers into a vector
	allcols = unlist(c(colnumbs, unlist(args[numb.ones])))

	#### return data without the columns
        if (keep){
            data = data[,allcols]
        } else {
            data = data[,-allcols]
	}
        
	return(data)
}
