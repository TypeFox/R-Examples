#' Reads the met.dat input file
#' 
#' @description Reads the meteorological input data in the met.dat file.
#' 
#' 
#' @param filename Default name of the met.dat file.
#' @param nlines Optional, how many lines of the metfile to read?
#' @return Returns a dataframe.
#' @author Remko Duursma
#' @keywords utilities
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' # Simple as pi:
#' metdata <- readmet()
#' 
#' }
#' 
#' @export
readmet <- function(filename="met.dat", nlines=-1){

	metlines <- readLines(filename)
	namesloc <- grep("^columns",tolower(metlines))
	namesline <- metlines[namesloc]
	datastart <- grep("DATA STARTS", metlines, ignore.case=TRUE)
	sp <- strsplit(namesline, "=")[[1]][2]
	NAMES <- delempty(str_trim(strsplit(sp, "\t")[[1]]))
	NAMES <- gsub("'", "", NAMES)
	NAMES <- unlist(strsplit(NAMES, " ", fixed=TRUE))
	metdata <- read.table(filename, header=FALSE, skip=datastart, nrows=nlines)
	names(metdata) <- NAMES
	names(metdata)[names(metdata) == "RH%"] <- "RHperc"

return(metdata)
}

