##' Extract only part of a string, given a separator
##'
##' Given a string with a separator (e.g., "Subject 001", where the separator is a space), this function can be used to extract only whatever follows the separator (in this case, "001").
##' It is often used when data comes in with a conglomorated identifier (such as case-matchNumber-drawNumber-Month).
##' @title Extract only part of a string
##' @param string a vector of strings
##' @param sep the separator that separates the parts of the strings
##' @param position the position of the element you wish to extract
##' @return the element the user wishes to extract
##' @author Dustin Fife
##' @export
##' @examples
##' barcode = c("Case-001-B", "Control-001-A", "Case-002-A")
##' subsetString(barcode, sep="-", position=2)
##' subsetString(barcode, sep="-", position=3)
subsetString = function(string, sep=" ", position=3){
	string = unlist(lapply(string, FUN=function(x){unlist(strsplit(x, sep, fixed=TRUE))[position]}))
	string	
}