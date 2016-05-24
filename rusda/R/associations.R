#' Downloads associations for input species from SMML Fungus-Host DB
#' 
#' Searches and downloads associations from SMML Fungus-Hosts Distributions and Specimens database
#' for fungus or plant species input vector
#' @param x a vector of class \code{character} containing fungal or plant species names or a genus name (see Details)
#' @param database a character string specifying the databases that should be queried. Valid are
#' \code{"FH"} (Fungus-Host Distributions), \code{"SP"} (Specimens) or \code{"both"} databases
#' @param spec_type a character string specifying the type of \code{x}. 
#' Can be either \code{"plant"} or \code{"fungus"}
#' @param clean logical, if \code{TRUE} a cleaning step is run of the resulting associations list
#' @param syn_include logical, if \code{TRUE} associations for synonyms are searched and added. For a
#' complete synonyms list check \code{rusda::synonyms}
#' @param process logical, if \code{TRUE} downloading and extraction process is displayed
#' 
#' 
#' @details The Fungus-Hosts distributions database 'FH' comprises data compiled from Literature. In
#' the uncleaned output all kinds of unspecified substrates are documented like "submerged wood".
#' Cleanded data displayes Linnean names only and species names with either "subsp.","f. sp." "f.",
#' "var.". The Specimens database comprises entries from field collections.
#' 
#' If genera names are supplied, then species are derived from the NCBI taxonomy.
#' 
#' 
#' @return an object of class \code{list}. 
#' @return First is synonyms, second is associations. Synonmys is a
#' vector of mode \code{list} with synonyms for \code{x}. Notice: This is not a
#' complete list of synonym data in the database. This is the list of synonyms that contain data for
#' the input \code{x}. For a complete synonyms list check \code{rusda::synonyms} or (if needed) for fungi R package rmycobank.
#' @return Associations is a vector of mode \code{list} of associations for \code{x}
#' 
#' @author Franz-Sebastian Krah
#' 
#' @examples
#' \dontrun{
#' ## Example for species name(s) as input
#' x <- "Fagus sylvatica"
#' pathogens <- associations(x, database = "both", clean = TRUE, syn_include = TRUE,
#' spec_type = "plant", process = TRUE)
#' x <- "Rosellinia ligniaria"
#' hosts <- associations(x, database = "both", clean = TRUE, syn_include = TRUE, 
#' spec_type = "fungus", process = TRUE)
#' is.element("Rosellinia ligniaria", pathogens$association[[1]])
#' is.element("Fagus sylvatica", hosts$association[[1]])
#' 
#' ## Example for genus/genera name(s) as input
#' x <- "Zehneria"
#' # or
#' x <- c("Zehneria", "Momordica")
#' hosts <- associations(x, database = "both", clean = TRUE, syn_include = TRUE, 
#' spec_type = "plant", process = TRUE)
#' }

associations <- function(x, database = c("FH", "SP", "both"), 
  spec_type = c("plant", "fungus"), clean = TRUE, syn_include = TRUE, process = TRUE)
{
  # test internet conectivity
  if(!url.exists("r-project.org") == TRUE) stop( "Not connected to the internet. Please create a stable connection and try again." )
  if(!is.character(getURL("http://nt.ars-grin.gov/fungaldatabases/index.cfm"))) stop(" Database is not available : http://nt.ars-grin.gov/fungaldatabases/index.cfm")
  # test if arguments are given
  expect_match(spec_type, ("fungus|plant"))
  expect_match(database, ("FH|SP|both"))
  
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
  
  ## I. PARSE DATA    ##
  ######################
  if(process == TRUE) { message("... retrieving data ... for:") }
  p <- foreach(i = seq_along(tax)) %do% getHF(tax[[i]], process, spec_type = spec_type)
  
  ## II. DATA CONDITIONS ##
  #########################
  taxa <- lapply(tax, function(x) { paste(as.character(x[1]), as.character(x[2])) })
  co <- lapply(p, getCOND)
  
  ## III. SYNONYMS ##
  ###################
  if(process == TRUE) { message("... extracting Synonyms ...") }
  syns <- lapply(p, getSYNS, process = process, taxa = taxa) 
  names(syns) <- taxa
  
  ## IV. EXTRACRING DATA  ##
  ##########################
  # FH DB #
  if(process == TRUE & database == "FH" | database == "both") { message("... extracting Fungus-Hosts DB ...") }
  i <- NULL
  hosts_hf <- foreach(i = seq_along(taxa)) %do%  {
    if(length(co[[i]]$hfu) == 0 | length(co[[i]]$hf.st) == 0){ hf <- "nodata" }
    if(length(co[[i]]$hf.st) > 0)
    {# Stop 
      hf.c <- grep("The Literature database has", p[[i]])
      ifelse(length(hf.c) > 0, hf.sp <- hf.c, 
        hf.sp <- (grep("No records were found in the Literature database", p[[i]])))
      if(length(hf.sp) == 0){
        hf.sp <- (grep(paste("There are no records for ",taxa[[i]], 
          " in the Literature database", sep=""), p[[i]]))}
      # extract
      p[[i]][(co[[i]]$hf.st + 1):(hf.sp - 1)] 
    }
  }
  names(hosts_hf) <- unlist(taxa)
  
  # Speciments DB  #
  if(process == TRUE & database == "SP" | database == "both") { message("... extracting Specimens DB ...") }
  i <- NULL
  hosts_sp <- foreach(i=seq_along(taxa)) %do% {
    if(length(co[[i]]$sp) == 0 | length(co[[i]]$spe.st) == 0){ specim <- "nodata" }
    if(length(co[[i]]$spe.st) > 0)
    {
      spe.sp <- grep("Systematic Mycology and Microbiology Laboratory[.]",p[[i]])
      spe.sp <- spe.sp[length(spe.sp)]
      specim <- p[[i]][(co[[i]]$spe.st + 1):(spe.sp - 1)]
    }
  }   
  names(hosts_sp) <- unlist(taxa)
  
  ## IV. SYNONYMS EXCLUDE  ##
  ###########################
  ## Exclude results for synonyms if wanted:
  # find occurences for taxon for stop condition and 
  # extract until next synonym of input taxa
  if(syn_include == FALSE){
    if(process == TRUE) { message("... excluding synonyms ...") }
    no_syns <- function(x){
      # search start and stops
      st <- foreach(i = seq_along(taxa)) %do% grep(taxa[[i]], x[[i]])
      sp <- foreach(i = seq_along(taxa)) %do% {
        sy <- paste(syns[[i]][!syns[[i]] == taxa[[i]]], collapse = "|")
        grep(sy, x[[i]], value = FALSE)}
      # choose next stop if there
      spp <- list(); for(i in seq_along(taxa)){
        if(is.integer(sp[[i]]) && length(sp[[i]]) == 0L){spp[[i]] <- integer(0)}
        if(length(st[[i]]) > 0 & length(sp[[i]]) > 0)    
          # choose the value next higher from starting (st) point, so two conditions 
          # must be matched: bigger and next integer, so the one with min distance
          spp[[i]] <- sp[[i]][(sp[[i]] > st[[i]][1]) & sp[[i]] == ((min(st[[i]][1] - sp[[i]]) * -1) + st[[i]][1])]
      }
      # choose only  
      res <- list(); for(i in seq_along(taxa)){
        # if there is no start and no stop
        if(length(st[[i]]) == 0 & is.integer(spp[[i]]) && length(spp[[i]]) == 0L){res[[i]] <- x}
        # if start but no stop: from stop to end (happens if input species occures at the bottom)
        if(length(st[[i]]) > 0 & length(spp[[i]]) == 0L)
          res[[i]] <- x[[i]][st[[i]][1] : length(x[[i]])]
        # if there is a start and stop condition
        if(length(st[[i]]) > 0 & length(spp[[i]]) > 0)
          if(length(st[[i]]) > 0)res[[i]] <- x[[i]][st[[i]][1]:(spp[[i]] - 1)]
          else
            res[[i]] <- x[[i]][st[[i]]:(spp[[i]] - 1)] }
      return(res)
    }
    res <- lapply(list(hosts_hf, hosts_sp), no_syns)
    hosts_hf <- res[[1]]
    hosts_sp <- res[[2]]
  }
  
  ## V. RESULTS OBJECT ##
  #######################
  if (database == "FH") { res <-  hosts_hf }
  if (database == "SP") { res <-  hosts_sp }
  if (database == "both") { 
    res <-  foreach(i = seq_along(hosts_hf)) %do% c(hosts_hf[[i]], hosts_sp[[i]])
    names(res) <- names(hosts_hf)
    res <- lapply(res, function(x)
    { if(length(grep("nodata",x)) == 2) { x <- "nodata" }
      if(!length(grep("nodata",x)) == 2) { x }})
  }
  
  ## VI. CLEAN    ##
  ##################
  ## do not conduct clean step if wanted
  if(clean == TRUE)
  {
    if(process == TRUE) { message("... cleaning step ...")}
    res <- lapply(res, clean_step, taxa = taxa, 
      syns = syns, spec_type = spec_type, synonyms_incl = TRUE)
  }
  res <- lapply(res, unique)
  names(res) <- taxa
  
  # VII. RESULTS OBJECT 2  ##
  ###########################
  return(list(synonyms = syns, associations = res))
}
