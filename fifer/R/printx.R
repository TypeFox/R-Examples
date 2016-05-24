##' Change Defaults for print.xtable
##'
##' @title Change Defaults of print.xtable
##' @details This is a wrapper for the function \code{\link{print.xtable}}, where the defaults have been specified
##' for apa-like tables. 
##' @param x The \code{xtable} object the user wishes to export
##' @param file the location where the file is stored
##' @param ... other arguments passed to print.xtable
##' @author Dustin Fife
##' @export
##' @import xtable
printx = function(x, file="",...){
	print(x, file=file, type="latex", 
		caption.placement="top", latex.environment="center", ...)		
}
