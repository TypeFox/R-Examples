#' Reads the dayflx.dat output file
#' 
#' @description Reads the dayflx.dat MAESTRA/MAESPA output file, returns a clean dataframe.
#' Names of the variables are read from the Columns: line.
#' 
#' 
#' @param filename Default name of the daily flux file.
#' @return Returns a dataframe.
#' @author Remko Duursma
#' @keywords utilities
#' @examples
#' 
#' 
#' \dontrun{
#' # Read it:
#' mysim1 <- readdayflux()
#' }
#' 
#' @export
readdayflux <- function(filename="dayflx.dat"){

	daylines <- readLines(filename, 100)
	colloc <- grep("Columns",daylines)
	dayflux <- read.table(filename, skip=colloc)
	names(dayflux) <- delempty(strsplit(delempty(strsplit(daylines[colloc], "Columns:")[[1]])," ")[[1]])
	 
	#names(dayflux) <- c("DOY","Tree","absPAR","absNIR","absTherm","totPs","totRf","netPs","totLE1","totLE2","totH")

dayflux
}

