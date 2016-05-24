#' Downloads and evaluate species presence in SMML DBs
#' 
#' Searches, downloads and evaluates presence/absence of data in the SMML databases
#' @param x a vector of class \code{character} containing fungal or plant species or genus names
#' @param spec_type a character string specifying the type of \code{x}. 
#' Can be either \code{"plant"} or \code{"fungus"}
#' @param process logical, if \code{TRUE} downloading and extraction process is displayed
#' 
#' @details Use this function before deriving data from one of the databases in order to prune your
#' input species vector. With pruned species vectors the functions will run faster. This is important
#' if \code{x} is some hundred species long. 
#' 
#' @return an object of class \code{data.frame}: presence/absence
#' 
#' @author Franz-Sebastian Krah
#' 
#' @examples
#' \dontrun{
#' fungus.meta <- meta_smml(x = "Picea abies", process = TRUE, spec_type = "plant")
#' fungus.meta
#' hosts.meta <- meta_smml(x = "Antrodiella citrinella", process = TRUE, spec_type = "fungus")
#' hosts.meta
#' }

meta_smml <- function(x, spec_type = c("plant", "fungus"), process = TRUE)
{
  if(!url.exists("r-project.org") == TRUE) stop( "Not connected to the internet. Please create a stable connection and try again." )
  if(!is.character(getURL("http://nt.ars-grin.gov/fungaldatabases/index.cfm"))) stop(" Database is not available : http://nt.ars-grin.gov/fungaldatabases/index.cfm")
  expect_match(spec_type, ("fungus|plant"))
  if(length(grep("\\sx\\s", x)) > 0) stop(" no hybrids allowed as input ")
  
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
  # produce taxon set
  if(length(grep(" ", x)) >=0) {tax <- strsplit(x, " ")}
  
  ## I. PARSE DATA       ##
  #########################
  if(process == TRUE) { message("... retrieving data ... for:") }
  p <- foreach(i = seq_along(tax)) %do% getHF(tax[[i]], process, spec_type = spec_type)
  i <- NULL
  
  # II. META Anylysis   ##
  ########################
  l <- lapply(p, getMETA)
  
  ## III. Create Object ##
  ########################
  names(l) <- x
  obj = do.call(rbind,l)
  return(obj)
}
