#' Reads the hrflux.dat MAESTRA/MAESPA output file
#' 
#' @description Reads the hourly output file (hrflux.dat).
#' 
#' 
#' @param filename Default name of the (half-)hourly output file.
#' @return Returns a dataframe.
#' @author Remko Duursma
#' @keywords utilities
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' # Simple as this:
#' mysim2 <- readhrflux()
#' }
#' @export
readhrflux <- function(filename="hrflux.dat"){
	
	hrlines <- readLines(filename, 100)
	colloc <- grep("Columns",hrlines)

	hrflux <- read.table(filename, skip=colloc,na.strings='NaN')
	names(hrflux) <- delempty(strsplit(delempty(strsplit(hrlines[colloc], "Columns:")[[1]])," ")[[1]])

	hrflux$conttime <- with(hrflux, DOY + HOUR/max(HOUR))
	
	hrflux
}

