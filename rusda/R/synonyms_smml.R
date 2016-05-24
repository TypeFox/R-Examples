#' Downloads synonym data from SMML Nomenclature DB
#' 
#' Searches and downloads synonym data from SMML Nomenclature database
#' @param x a vector of class \code{character} containing fungal or plant species or genus names
#' @param spec_type a character string specifying the type of \code{x}. 
#' Can be either \code{"plant"} or \code{"fungus"}
#' @param clean logical, if \code{TRUE} a cleaning step is run of the resulting associations list
#' @param process logical, if \code{TRUE} downloading and extraction process is displayed
#' 
#' @return an object of class \code{list} containing synonyms for \code{x}
#' 
#' @author Franz-Sebastian Krah
#' 
#' @examples
#' \dontrun{
#' x <- "Solanum tuberosum"
#' synonyms_usda(x, spec_type = "plant", process = TRUE, clean = TRUE)
#' x <- c("Phytophthora infestans", "Polyporus badius")
#' synonyms_usda(x, spec_type = "fungus", process = TRUE, clean = TRUE)
#' }

synonyms_smml <- function(x, spec_type = c("plant", "fungus"), clean = TRUE, process = TRUE)
{
  if(!url.exists("r-project.org") == TRUE) stop( "Not connected to the internet. Please create a stable connection and try again." )
  if(!is.character(getURL("http://nt.ars-grin.gov/fungaldatabases/index.cfm"))) stop(" Database is not available : http://nt.ars-grin.gov/fungaldatabases/index.cfm")
  expect_match(spec_type, ("fungus|plant"))
  if(length(grep("\\sx\\s", x)) > 0) { stop(" no hybrids allowed as input") }
  
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

  if(length(grep(" ", x)) >=0) {tax <- strsplit(x, " ")}
  i <- NULL
  
  ## I. PARSE DATA    ##
  ######################
  if(process == TRUE) { message("... retrieving data ... for:") }
  p <- foreach(i = seq_along(tax)) %do% getHF(tax[[i]], process = process, spec_type = spec_type)
  taxa <- lapply(tax, function(x) { paste(as.character(x[1]), as.character(x[2])) })
  
  ## III. SYNONYMS ##
  ###################
  if(process == TRUE) { message("... extracting Synonyms ...") }
  syns <- lapply(p, getSYNS, process = process, taxa = taxa)
  synos <- list(); for(i in seq_along(taxa))
  {
    st <- grep(paste("Nomenclature data for", as.character(taxa[[i]]), sep = " "), p[[i]])
    cond <- paste(c("Distribution", "Substrate", "Updated", "The Fungus-Host Distributions", "Notes"), collapse="|")
    sp <- grep(cond, p[[i]])
    if(length(st) == 0){ synos[[i]] <- "nodata" }
    else{ synos[[i]] <- p[[i]][(st + 1):(sp[1] - 1)] }
  }
  synos <- foreach(i = seq_along(taxa)) %do% c(syns[[i]], synos[[i]])
  
  if(clean == TRUE)
  {
    synos <- lapply(synos, clean_step, taxa = taxa, spec_type = spec_type, syns = syns, synonyms_incl = FALSE)
  }
  synos <- lapply(synos, unique)
  names(synos) <- taxa 
  return(synos)
}
