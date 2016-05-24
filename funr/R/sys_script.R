


## http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script#comment12780420_1815606

#' Get path of current script
#'
#' @export
#'
#' @source
#' http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
#'
#' @examples \dontrun{
#' ## cat foo.R
#' ## library(funr);sys.script()
#'
#' ## source("foo.R")
#'
#' ## Rscript foo.R
#'
#' }
sys.script <- function(){
	cmdargs <- commandArgs(trailingOnly = FALSE)
	fl <- grep("--file=", cmdargs)
	if (length(fl) > 0) {
		# Rscript
		return(normalizePath(gsub("--file=", "", cmdargs[fl])))
	} else {
		# 'source'd via R console
		return(normalizePath(sys.frames()[[1]]$ofile))
	}
}
