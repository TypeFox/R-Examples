#' Replaces a namelist or a parameter
#' 
#' @description The function \code{replaceNameList} replaces the whole namelist in an input file. 
#' All parameters in the namelist must be provided, otherwise MAESTRA/MAESPA will likely crash. Or, you can use
#' \code{replacePAR} to replace a single parameter. If the new parameter value(s) is a vector (or a single value), the values #' will be placed on a single line in the parameter file. If instead a matrix is provided, each row of the matrix is placed on a separate line. 
#' 
#' WARNING : The input file is modified. Make sure to backup your original
#' input files!
#' 
#' 
#' @param namelist Name of the namelist.
#' @param datfile Name of the input file.
#' @param vals A list of values (see example below).
#' @param parname Name of the parameter to replace the value of.
#' @param newval New value of the parameter. Can be a single value or a vector,
#' or a matrix (see Details).
#' @param noquotes Logical. If FALSE, does print quotes around character
#' values.
#' @return Nothing is returned. The input file is modified.
#' @author Remko Duursma
#' @seealso \code{\link{replacePAR}}
#' @keywords utilities
#' @examples
#' 
#' 
#' \dontrun{
#' # Replace an entire namelist
#' replaceNameList(namelist="aerodyn", 
#'     datfile="trees.dat", vals=list(zht=30,zpd=3,z0ht=0.6))
#' 
#' #' # Replace a parameter with a single value:
#' replacePAR("trees.dat", "notrees", "plot", newval=100)
#' 
#' # Replace a number of values:
#' replacePAR("trees.dat", "xycoords", "xy", newval=c(1,1,2,2,3,3))
#' 
#' # Replace a number of values in a different way : this may be useful in 
#' # more complicated programs,
#' # rr when reading a string from a file (that should be interpreted as a vector).
#' replacePAR("trees.dat", "xycoords", "xy", newval="1 1 2 2 3 3", noquotes=TRUE)
#' 
#' # Replace a parameter so that the new values are placed on multiple rows.
#' # Useful for tree namelists with multiple dates and multiple trees 
#' # (where each row corresponds to a tree, and each column to a date.)
#' m <- matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, byrow=TRUE)
#' replacePAR("trees.dat", "values", "indivlarea", newval=m)
#' }
#' 
#' @export
#' @rdname replacePAR
replaceNameList <- function(namelist, datfile, vals){
    
    # Find NAMELIST and end of it ("/")
    datlines <- str_trim(readLines(datfile))
    nmreg <- paste0("&",namelist,"$")
	
    nl_start <- grep(nmreg, datlines, ignore.case=TRUE)
    
    # Find nearest namelist closer ("/")
    endnml <- grep("^/$", datlines)
    nmllen <- min(endnml[endnml > nl_start] - nl_start)
    nl_end <- nl_start + nmllen
    
    datlines_namelist <- datlines[nl_start:(nl_start + nmllen)]
    
    # Part of original file before and after this namelist.
    if(nl_start > 1){
      pref <- datlines[1:(nl_start-1)]
    } else {
      pref <- ""
    }
    if((nl_end+1) < length(datlines)){
      postf <- datlines[(nl_end+1):length(datlines)]
    } else {
      postf <- ""
    }
    
    # New namelist
    newlist <- formatNameList(namelist, vals)    
    
    # Rewrite file
    Lines <- c(pref, newlist, postf)
    writeLines(Lines, datfile) 
}

 
formatNameList <- function(namelist, vals){
  
  # New namelist
  listvals <- c()
  for(i in 1:length(vals)){
    
    if(!is.matrix(vals[[i]])){
      listvals[i] <- paste(names(vals)[i],printme(vals[[i]]), 
                         sep=" = ")
    } else {
      m <- vals[[i]]
      listvals[i] <- paste(names(vals)[i], paste(paste(apply(m,1,paste,collapse=" "), collapse="\n"),"\n"),
                           sep= " = ")
    }
    newlist <- c(paste0("&",namelist),
                 listvals, "/")
  }

return(newlist)
}




#'@rdname replacePAR
#'@export
replacePAR <- function(datfile, parname, namelist, newval, noquotes=TRUE){

  p <- readNameList(datfile, namelist)
  
  # Quoting or not? Do not listen to argument when newval is a matrix.
  if(!noquotes && !is.matrix(newval))newval <- printme(newval)
  
  p[[parname]] <- newval
  
  replaceNameList(namelist, datfile, p)
}





