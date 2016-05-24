#'set run options for running ADMB via R
#'
#'A helper function
#'
#'
#'@param check_tpl Check the specified TPL file for problems?
#'@param write_files Write out data and initialization files?
#'@param checkparam How to check \code{PARAMETERS} section of the TPL file:
#'\code{stop}=stop if there are problems; \code{warn}=give a warning if there
#'are problems, but try to proceed; \code{write}=modify TPL file, writing
#'appropriate sections; \code{ignore}=assume TPL file is OK, proceed
#'@param checkdata as with \code{checkparam}: how/whether to check/generate the
#'\code{DATA} section of the TPL file
#'@param compile compile the TPL file (via ADMB) into an executable?
#'@param run run the executable file with the specified data/initial values?
#'@param read_files read the results of an ADMB run into R?
#'@param clean_files Delete working files after completion of the run?  Options
#'are \code{"all"}, \code{"sys"}, \code{"output"}, \code{"none"}; \code{TRUE}
#'is equivalent to \code{"all"} and \code{FALSE} is equivalent to \code{"none"}
#'
#'@return A list with appropriate default values inserted for passing to
#'\code{\link{do_admb}}
#'@export
#'@author Ben Bolker
#'@keywords misc
run.control <- function(check_tpl=TRUE,
		write_files=TRUE,
		checkparam=c("stop","warn","write","ignore"),
		checkdata=c("stop","warn","write","ignore"),
		compile=TRUE,
		run=TRUE,
		read_files=TRUE,
		clean_files="all") {
	checkparam <- match.arg(checkparam)
	checkdata <- match.arg(checkdata)
	if (is.logical(clean_files)) {
		if (clean_files) clean_files <- "all"
		else clean_files <- "none"
	}
	c(check_tpl=check_tpl,checkparam=checkparam,checkdata=checkdata,
			write_files=write_files,compile=compile,run=run,read_files=read_files,
			clean_files=clean_files)
}

