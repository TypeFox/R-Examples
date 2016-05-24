#' Downloads literature from SMML Literature DB
#' 
#' Searches and downloads literature entries from the SMML Literature database
#' 
#' @param x a vector of class \code{character} containing fungal or plant species names
#' @param spec_type a character string specifying the type of \code{spec}. Can be either 
#' \code{"plant"} or \code{"fungus"}
#' @param process logical, if \code{TRUE} downloading and extraction process is displayed
#' 
#' an object of class \code{list} 
#' @return a vector of mode \code{list} with literature entries for \code{x}
#' 
#' @author Franz-Sebastian Krah
#' 
#' @examples
#' \dontrun{
#' x <- "Polyporus badius"
#' lit <- literature(x, process = TRUE, spec_type = "fungus")
#' lit
#' }

literature <- function(x, spec_type = c("plant", "fungus"), process = TRUE)
{
  if(!url.exists("r-project.org") == TRUE) stop( "Not connected to the internet. Please create a stable connection and try again." )
  if(!is.character(getURL("http://nt.ars-grin.gov/fungaldatabases/index.cfm"))) stop(" Database is not available : http://nt.ars-grin.gov/fungaldatabases/index.cfm")
  expect_match(spec_type, ("fungus|plant"))
  if(length(grep("\\sx\\s", x)) > 0){ stop(" no hybrids allowed ") }
  if(length(spec_type) == 2) stop(" 'spec_type' not specified. Please choose one of 'plant', 'fungus'")
  ifelse(length(grep(" ", x)) > 0,tax <- strsplit(x, " "), tax <- strsplit(x, "_"))
  
  ## I. PARSE DATA    ##
  ######################
  if(process == TRUE) { message("... retrieving data ... for:") }
  p <- foreach(i = seq_along(tax)) %do% getHF(tax[[i]], process = process, spec_type = spec_type)
  
  ## II. DATA CONDITIONS ##
  ######################### 
  taxa <- lapply(tax, function(x){paste(as.character(x[1]), as.character(x[2]))})
  co <- lapply(p, getCOND)
  
  ## IV. DATA EXTRACTION ##
  #########################
  i <- NULL
  l <- foreach(i = seq_along(p)) %do% 
  {
    l.st <- grep("The Literature database has" , p[[i]])
    # Stop
    ifelse(length(l.sp <- grep("The Specimens database has", p[[i]])) > 0,l.sp,
    l.sp <- grep(paste("There are no records for",taxa[[i]], "in the Specimens database"), p[[i]]))
    
    lit <- p[[i]][(l.st + 1):(l.sp - 1)] 
    if(length(which(nchar(lit)==0)) > 0) {lit[-which(nchar(lit) == 0)] }
    else lit
  }
  return(l)
}
