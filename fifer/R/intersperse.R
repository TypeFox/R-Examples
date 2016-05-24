##' Intersperse elements of two+ vectors
##' Given two vectors (one of length i, the other of length j), \code{intersperse} will combine the elements of each vector into strings of length
##' i X j, where each element is the concatenation of the elements of the two vectors. See examples.  
##' @title Intersperse elements of Two+ Vectors
##' @param ... the vectors the user wishes to intersperse
##' @return a vector of length i X j, containing the interspersed vectors as strings
##' @author Dustin Fife
##' @export
##' @examples
##' intersperse(LETTERS[1:3], 1:3)
intersperse = function(...){
	nn = expand.grid(..., stringsAsFactors=F)
	do.call("paste0", c(nn, sep=""))
}