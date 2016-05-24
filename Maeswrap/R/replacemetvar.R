#' Replace a weather variable
#' 
#' Replaces one (or more) of the weather variables in the met.dat file.
#' 
#' 
#' @aliases replacemetvar replacemetdata
#' @param replacevar Character. Name(s) of the variable to be replaced.
#' @param newvalues Vector of new values for the weather variable, has to be
#' the same length as the number of records in the met.dat file.
#' @param oldmetfile Default name of the met.dat file that will be modified.
#' @param newmetfile Name of the new met.dat file.
#' @param metdfr Dataframe with met data, to be pasted into a met.dat file.
#' @param columns Optional character string : if the 'Columns' statement in the
#' met.dat file is to be replaced.
#' @param khrs Optional. Number of timesteps per day (by default, read from the
#' met.dat file).
#' @param setdates Logical. If TRUE, fixes the start and end date in the new met file.
#' @return Returns nothing.
#' @author Remko Duursma
#' @keywords utilities
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' #:::1.::: Replace precipitation with random number between 0 and 2.
#' # First find out how many records there are:
#' nrecords <- nrow(readmet("met.dat"))
#' 
#' # Make new rain
#' newrain <- runif(nrecords, 0, 2)
#' 
#' # And replace
#' replacemetvar("PPT",newrain,"met.dat", "newmet.dat")
#' 
#' 
#' #:::2.::: Replace multiple weather variables.
#' newtair <- runif(nrecords, 0, 35)
#' 
#' # Have to make a matrix of the variables to be replaced:
#' newmat <- matrix(cbind(newrain, newtair),ncol=2)
#' 
#' # And give a vector of variable names --in the same order as in the matrix!!--.
#' replacemetvar(c("PPT","TAIR"), newmat, "met.dat", "newmet.dat")
#' 
#' 
#' }
#' 
#' @rdname replacemetvar
#' @importFrom utils write.table
#' @export
replacemetvar <- function(replacevar, newvalues, oldmetfile="met.dat", newmetfile="metNEW.dat"){
	metlines <- readLines(oldmetfile)
	datastart <- grep("DATA START", metlines, ignore.case=TRUE)
	DATA <- read.table(oldmetfile, skip=datastart, header=FALSE)    
	preamble <- readLines(oldmetfile)[1:datastart]
  namesloc <- grep("columns", metlines, ignore.case=TRUE)
  namesloc <- setdiff(namesloc, grep("nocolumns", metlines, ignore.case=TRUE))
  namesline <- metlines[namesloc]
  sp <- strsplit(namesline, "=")[[1]][2]
  NAMES <- delempty(strsplit(sp, "\t")[[1]])
  NAMES <- str_trim(gsub("'", "", NAMES))
  NAMES <- do.call("c", strsplit(NAMES, " ", fixed = TRUE))
  names(DATA) <- NAMES
    
	if(nrow(DATA) != length(newvalues) | (is.matrix(DATA) && nrow(DATA) != nrow(newvalues)))
        stop("Length of new data not equal to nr rows in metfile")

  DATA[,replacevar] <- newvalues
	
	writeLines(preamble,newmetfile)
	write.table(DATA, newmetfile, sep=" ",row.names=FALSE,col.names=FALSE,append=TRUE)
} 
