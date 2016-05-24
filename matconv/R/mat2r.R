#' mat2r
#'
#'The top level driver function to call the converting functions and handle the
#'input and output.
#' @param inMat A file path with the input Matlab / Octave code to be converted
#' or a character vector of the code that needs to be converted
#' @param pathOutR A file path with the desired output file location
#' @param funcConverters A list of function converters that the user wants to
#'   use in this conversion made by \code{\link{makeFuncMaps}}
#' @param dataConverters A list of data converters that the user wants to 
#'   use in this conversion made by \code{\link{makeSliceMap}} or
#'   \code{\link{makeDataMap}}
#' @param verbose A number indicating the amount of messages that should be outputed.
#' \describe{
#'   \item{0}{No messages}
#'   \item{1}{A summary report of what happened in the conversion}
#'   \item{2}{The final code as a message as well as the summary report}
#' }
#'
#' @return A list containing the original code (named matCode) and the converted code (named rCode).
#' @examples
#' matIn <- c("function [ dat ] = xlsReadPretty(varargin)", 
#'  "  didThing = 1*3;",
#'  "  dat = didThing / 3;",
#'  "end")
#'  mat2r(matIn, verbose = 0)$rCode
#' 
#' # "xlsReadPretty <- function(...){" 
#' # "\tvarargin <- list(...)"
#' # "  didThing <- 1*3"
#' # "  dat <- didThing / 3"
#' #"\treturn(dat)"
#' #"}"
#' @export
mat2r <- function(inMat,
                  pathOutR ='',
                  funcConverters = NULL,
                  dataConverters = NULL,
                  verbose = 1){

	if (length(inMat) == 1 && file.exists(inMat)){
		if(!grepl("[.]m", inMat)) stop("Please supply a '.m' file")
		rawMat <- readLines(inMat)
	} else {
		rawMat <- inMat
	}

	linesDes <- linesOrg <- trimWhite(rawMat, "end")
	isScr <- !grepl("function", linesOrg[1])



	commentSet <- grepl("^%", linesDes) | grepl("^\\s+%", linesDes)

	codeDes <- linesDes[!commentSet]
	codeDes <- convEasySyntax(codeDes)

	linesDes[!commentSet] <- codeDes
	linesDes[commentSet] <- gsub("%", "#", linesDes[commentSet])

	if (!isScr) linesDes <- convUserFunctions(linesDes)

	if(!is.null(dataConverters)){
		linesDes <- convData(linesDes, dataConverters)
	}

	if(!is.null(funcConverters)){
		linesDes <- convFunctionsCalls(linesDes, funcConverters)
	}



	report <- sprintf("The previous code had %d lines and the R code has %d lines",
                    length(linesOrg),
                    length(linesDes))

	if(verbose == 2 ){
		message(report)
		message(linesDes)
	} else if (verbose == 1) {
		message(report)
	}

	if(nzchar(pathOutR)) writeLines(linesDes, pathOutR)

	return(list(matCode = linesOrg, rCode = linesDes))

}
