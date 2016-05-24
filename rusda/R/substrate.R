#' Downloads substrate data from SMML Nomenclature DB
#' 
#' Searches and downloads substrate data from SMML Nomenclature database
#' 
#' @param x a vector of class \code{character} containing fungal or plant species names
#' @param process logical, if \code{TRUE} downloading and extraction process is displayed
#' 
#' @details Don't be disappointed. Not much data there. 
#' But depends on the study group, so give it try.
#' 
#' @return an object of mode \code{list} containing substrate for fungus species
#' 
#' @author Franz-Sebastian Krah
#' 
#' @examples
#' \dontrun{
#' x <- c("Polyporus_rhizophilus", "Polyporus_squamosus")
#' subs.poly <- substrate(x, process=TRUE)
#' subs.poly
#' }
#' 

substrate <- function(x, process = TRUE)
{
  if(length(grep("\\sx\\s", x)) > 0) stop(" no hybrids allowed as input ")
  if(!url.exists("r-project.org") == TRUE) stop( "Not connected to the internet. Please create a stable connection and try again." )
  if(!is.character(getURL("http://nt.ars-grin.gov/fungaldatabases/index.cfm"))) stop(" Database is not available : http://nt.ars-grin.gov/fungaldatabases/index.cfm")
  
  # if underscore remove it
  if(length(grep("_", x)) > 0) {x <- gsub("_", " ", x)}
  
  ## If a genus is supplied
  words <- vapply(strsplit(x, "\\W+"), length, integer(1))
  if (any(words == 1) & any(words == 2)) 
  {stop(paste(" check if you specified ONLY genus names or ONLY species names \n",
    "AFAICS you provided:  \n", sum(words==1), "  genus name(s)  \n", sum(words==2), "  species name(s) ", sep=""))}
  if(all(words == 1)){
    x <- lapply(x, ncbiSpecies, clean = TRUE, sub = FALSE)
    x <- unlist(x)
  }
  
  # tests
  if(length(grep("\\sx\\s", x)) > 0) 
    stop(" no hybrids allowed as input ")
  
  tax <- strsplit(x, " ")
  
  ## I. PARSE DATA       ##
  #########################
  if(process == TRUE) { message("... retrieving data ... for:") }
  p <- foreach(i = seq_along(tax)) %do% getHF(tax[[i]], process = process, spec_type = "fungus")
  i <- NULL
  
  ## II. SUBSTRATE       ##
  ######################### 
  subst <- function(x)
  {
    subs.st <- grep("Substrate:", x)
    subs <- x[subs.st]
    subs <- gsub("Substrate: ", "", subs)
    ifelse(length(subs.st) > 0, subs, "nodata")
  }
  obj <- lapply(p, subst)
  
  ## III. Create Object ##
  ########################
  names(obj) <- x
  return(obj)
}
