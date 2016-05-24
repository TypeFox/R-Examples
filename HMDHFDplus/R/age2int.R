
#'
#' @title \code{age2int()} convert the Age column from standard HMD or HFD tables to integer
#' 
#' @description Long the bane of many an HMD/HFD user is that the age column must be read into R as a factor or character vector, yet we'd like to use it as integer or numeric. This function strips symbols that are used to indicate the open age groups ("12-","55+","110+"), and coerces to integer format. This function is called by \code{HFDparse()} and \code{HMDparse()}, and so forth.
#' 
#' @param Age a vector of the Age column from and HMD or HFD data object that has been read directly into R. This may be a factor or character vector.
#' 
#' @export 
#' 
#' @return the same age vector as a clean integer.
#' 
#' @details This function is written for the sake of various parse functions.
#' 
#' @note original function submitted by Josh Goldstein, modified by Tim Riffe.
#' 
#' @examples 
#' AgeTest <- c("12-","13","14","55+")
#' (AgeNew  <- age2int(AgeTest))
#' AgeNew + .5 # sort of mid-interval
#' 
#' # also handles abrdiged ages properly:
#' AgeAbridged <- c("0","1-4","5-9","10-14")
#' age2int(AgeAbridged)
age2int <- function(Age){
    ## replaces + and - with nothing
    ## e.g., "99+" becomes "99"
    as.integer(gsub("[-+]", "", unlist(lapply(strsplit(as.character(Age),split="-"),"[[",1))))
}



