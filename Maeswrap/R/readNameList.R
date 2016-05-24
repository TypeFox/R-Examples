#' Read a namelist into a list
#' 
#' @param datfile Name of the Maestra/Maespa input file (e.g. "trees.dat")
#' @param namelist Name of the namelist (e.g. "plot")
#' @export
#' @rdname readNameList
readNameList <- function(datfile, namelist){
  
  r <- str_trim(readLines(datfile))
  
  nmStart <- grep(paste0("&",namelist, "$"), r, ignore.case=TRUE)
  r <- r[nmStart[1]:length(r)]
  r <- r[1:grep("^/$",r)[1]]
  r <- r[-c(1,length(r))]
  r <- delempty(r)
  
  # figure out which elements belong to which parameter.
  parloc <- grep("=",r)
  last <- length(r) - parloc[length(parloc)] + 1
  nlines <- c(diff(parloc),last)
  
  l <- list()
  for(i in 1:length(parloc)){
    ind <- parloc[i]:(parloc[i] + nlines[i] -1)
    l[[i]] <- parsePARline(r[ind])
  }
  
  parnames <- str_trim(sapply(strsplit(r[parloc],"="), "[", 1))
  names(l) <- tolower(parnames)
  
  return(l)
}

