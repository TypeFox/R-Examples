#' Reads the watbal.dat MAESPA output file
#' 
#' @description Reads the hourly water balance output file ("watbal.dat").
#' 
#' 
#' @param filename Default name of the (half-)hourly water balance output file.
#' @return Returns a dataframe.
#' @author Remko Duursma
#' @keywords utilities
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' # Simple as this:
#' mywatbalresult <- readwatbal()
#' 
#' # If you want to select the water balance file with a menu:
#' readwatbal(file.choose())
#' 
#' }
#' @export
#' @importFrom utils read.table
readwatbal <- function(filename="watbal.dat"){

	watlines <- readLines(filename, 100)
	colloc <- grep("Columns",watlines)
	namesline <- watlines[colloc]
	NAMES <- delempty(strsplit(strsplit(namesline, ":")[[1]][2], " ")[[1]])
	watbal <- read.table(filename, header=FALSE, na.strings="-999.0000", skip=colloc)
	names(watbal) <- NAMES
	
return(watbal)
}
