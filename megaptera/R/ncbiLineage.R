ncbiLineage <- function(taxon, species.list = FALSE, 
                        kingdom, parallel, retmax = 500){
  
  ## format 'taxon'
  ## --------------
  if ( species.list ){
    cat("\n.. taxon: ", paste(head(taxon, 2), collapse = ", "), 
        ", ... [", length(taxon), " species]", sep = "")
  } else {
    cat("\n.. taxon:", taxon, "..")
  }
  
  rank <- ifelse(species.list, 
                 "genus [rank]",
                 "species [rank] AND specified [prop]")
  taxon <- paste(taxon, "[subtree] AND", rank)
  if ( length(taxon) > 1 ) {
    taxon <- paste("(", taxon, ")", sep = "")
    taxon <- paste(taxon, collapse = " OR ")
  }
  
  ## get UID of 'taxon'
  ## ------------------
  cat("\n.. posting UIDs on Entrez History Server ..")
  xml <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/", 
               "eutils/esearch.fcgi?",
               "tool=megaptera",
               "&email=heibsta@gmx.net",
               "&usehistory=y", 
               "&db=taxonomy", 
               "&term=", taxon, sep = "")
  xml <- xmlParse(xml)
  
  ## get and parse results via eFetch
  ## --------------------------------
  
  webEnv <- xpathSApply(xml, fun = xmlToList,
                        path = "//eSearchResult/WebEnv")
  queryKey <- xpathSApply(xml, fun = xmlToList,
                          path = "//eSearchResult/QueryKey")
  
  ## get complete lineage for each species
  ## -------------------------------------
  cat("\n.. retrieve classification ..")
  retstart <- 0
  i <- 1
  x <- list()
  repeat {
    
    cat("\n   - processing batch", i)
    
    ## get XML with full records
    ## -------------------------
    #     slog("\n.. retrieving full records ( batch", i, ") ..", file = logfile)
    xml <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
                 "efetch.fcgi?tool=megaptera&email=heibsta@gmx.net",
                 "&db=taxonomy", "&query_key=", queryKey, "&WebEnv=", webEnv,
                 "&retmode=xml&retmax=", retmax, "&retstart=", retstart, sep = "")
    xml <- xmlParse(xml, getDTD = FALSE)
    newXMLNamespace(xmlRoot(xml), "http://www.w3.org/XML/1998/namespace", prefix = "ncbi", set = TRUE)
    #     xml <- xmlRoot(xml)
    #  xmlAttrs(xml) <- c(xmlns = "http://www.w3.org/XML/1998/namespace")
    
    ## check for exsitance of taxon name in both Metazoa and 
    ## Viridiplantae, e.g. Glaucidium, Prunella, etc.
    ## Currently, the check is done further downstream
    ## Here it makes sense only when disambiguation
    ## can be achieved via eSeach/Entrez History
    ## ----------------------------------------------
    #     kingdomCheck <- xpathSApply(xml, "//Taxon[Rank='kingdom']/ScientificName", xmlValue)
    #     kingdomCheck <- unique(kingdomCheck)
    
    spec <- xpathSApply(xml, "/ncbi:TaxaSet/Taxon/ScientificName", xmlValue)
    id <- grep("'", spec)
    dropped <- length(id)
    if ( dropped > 0 ){
      cat("      ignored unconventional species name:", spec[id])
      spec <- spec[-id] ## e.g., Gynnostemma 'burmanicum'
    }
    
    ##
    if ( parallel ) {
      # sfLibrary(megaptera)
      sfExport("xml")
      sfInit(parallel = parallel, cpus = 12, type = "SOCK")
      tmp <- sfLapply(spec, getLineage, xml = xml, kingdom = kingdom)
      sfStop()
    } else {
      tmp <- lapply(spec, getLineage, xml = xml, kingdom = kingdom)
    }
    ## kingdom ambiguity can cause NULL elements!
    x <- c(x, tmp[!sapply(tmp, is.null)])
    
    ## check if there are UIDs left to fetch
    ## -------------------------------------
    if ( length(tmp) == retmax - dropped){
      retstart <- length(x)
      i <- i + 1
    } else {
      break
    } 
  } # END OF WHILE-loop
  x
}