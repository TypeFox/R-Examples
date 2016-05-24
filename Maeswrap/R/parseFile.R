#' Parse an input file
#' 
#' @description Takes an input file for MAESTRA/MAESPA, and reads all namelists into a
#' nested list. Also reads the first line of the file, which (optionally) contains a title, 
#' to be used in Maestra/pa output files.
#' 
#' 
#' @param fn Filename
#' @return Returns a named list, each element contains a namelist and its
#' parameters.
#' @seealso To read one namelist from a file, see \code{\link{readNameList}}.
#' @examples
#' 
#' \dontrun{
#' # Parse a file
#' con <- parseFile("confile.dat")
#' 
#' # Namelists in the file
#' names(con)
#' }
#' 
#' @export parseFile
#' @importFrom stringr str_trim
parseFile <- function(fn){
  
  r <- str_trim(readLines(fn))
  
  Title <- r[1]
  
  nml <- r[grep("^&",r)]
  nml <- gsub("&","",nml)
  
  l <- list()
  for(i in 1:length(nml)){
    
    l[[i]] <- readNameList(fn, nml[i])
    
  }
  l <- c(Title,l)
  names(l) <- c("Title", nml)
  
return(l)
}
