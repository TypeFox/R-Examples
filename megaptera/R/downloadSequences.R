# DOWNLOAD SEQUENCES FROM GENBANK
# called by: stepAP
# package: megaptera
# author: Christoph Heibl (at gmx.net)
# last update: 2014-11-03

downloadSequences <- function(x, taxon){
  
  ## DEFINITIONS
  ## -----------
#   maxGIPerSpec <- x@maxGIPerSpec # currently not used
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  cols <- paste(gene, c("gb", "sel", "blocks"), sep = "_")
  logfile <- paste(acc.tab, "stepAP.log", sep = "-")
  retmax <- 500
  
  slog("\n.. search taxon:", taxon, "..", file = logfile)
  
  ## post UIDs on Entrez History Server
  ## ----------------------------------
  slog("\n.. posting UIDs on Entrez History Server ..", file = logfile)
  xml <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/", 
               "eutils/esearch.fcgi?",
               "tool=megaptera",
               "&email=heibsta@gmx.net",
               "&usehistory=y", 
               "&db=nucleotide", sep = "")
  xml <- paste(xml, "&term=", 
               term(taxon, x@locus),
               sep = "")
  
  ## get and parse results via eFetch
  ## --------------------------------
  xml <- xmlTreeParse(xml, useInternalNodes = TRUE)
  webEnv <- xpathSApply(xml, fun = xmlToList,
                        path = "//eSearchResult/WebEnv")
  queryKey <- xpathSApply(xml, fun = xmlToList,
                          path = "//eSearchResult/QueryKey")
  
  retstart <- 0
  i <- 1
  repeat {
    
    ## get XML with full records
    ## -------------------------
    slog("\n.. retrieving full records ( batch", i, ") ..", file = logfile)
    xml <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
                 "efetch.fcgi?tool=megaptera&email=heibsta@gmx.net",
                 "&db=nucleotide&query_key=", queryKey, 
                 "&WebEnv=", webEnv,
                 "&rettype=gb&retmode=xml",
                 "&retstart=", retstart, 
                 "&retmax=", retmax, sep = "")
    
    ## parse XML with xPath
    ## --------------------
    slog("\n.. parsing XML ..", file = logfile)
    xml <- xmlTreeParse(xml, getDTD = FALSE, useInternalNodes = TRUE)
    #     saveXML(xml, "megapteraAP-DEBUG.xml"); system("open -t megapteraAP-DEBUG.xml")
    ERROR <- xpathSApply(xml, fun = xmlValue, path = "//ERROR")
    if ( length(ERROR) > 0 ){
      if ( ERROR == "Empty result - nothing to do" ){
        slog("\n.. no sequences available ..", file = logfile)
        return(NULL)
      } else stop("unknown error in downloadSequences()")
    }
    
    ## get GI
    ## ------
    #     gi <- xpathSApply(xml, "//GBSeqid[matches(., 'gi')]", xmlValue)
    ## the preceding line should be preferred, but xPath2 doesn't seem
    ## to be implemented in libxml2 ???
    gi <- xpathApply(xml, "//GBSeq_other-seqids", xmlToList)
    gi <- lapply(gi, unlist)
    gi <- sapply(gi, function(x) x[grep("^gi[|]", x)])
    gi <- gsub("^gi[|]", "", gi)
    
    #     if ( any("1731829" %in% gi) ) stop(i)
    
    ## get organism = taxon name
    ## -------------------------
    organism <- xpathSApply(xml, "//GBSeq_organism", xmlValue)
    if ( length(grep(" " , organism)) == 0 )
      organism <- paste(organism, "sp.") # sometimes only a genus name is given in the organism field
    organism <- gsub(" ", "_", organism) # use underscores
    #     organism <- gsub("'|:", "", organism) # do not use special characters
    #     organism <- organism[1:length(dna)] # dirty hack!
    
    ## TOPOLOGY + LENGTH
    ## sequences > 200000 bp seem not to be included in XML!
    ## -----------------------------------------------------
    topology <- xpathSApply(xml, "//GBSeq_topology", xmlValue)
    bp <- as.numeric(xpathSApply(xml, "//GBSeq_length", xmlToList))
    contig <- length(xpathApply(xml, "//GBSeq_contig", xmlToList))
    ## dangerous hack!!! (Cetacea, RAG1)
    bp.id <- bp >= sort(bp, decreasing = TRUE)[contig]   ## 750000
    if ( any(bp.id) ){
      gi <- gi[!bp.id]
      organism <- organism[!bp.id]
      topology <- topology[!bp.id]
      bp <- bp[!bp.id]
    }  
    
    ## get DNA
    ## -------
    dna <- xpathSApply(xml, "//GBSeq_sequence", xmlToList)
    
    ## extracted locus from annotated genome
    ## -------------------------------------
    ## Note: We cannot rely on topology only,
    ## as sometimes people mark single cp genes
    ## as "circular" 
    id <- which( topology %in% "circular" & bp > 10000 )
    if ( length(id) > 0 ){
      slog(paste("\n.. extracting '", x@locus@aliases[1], 
                 "' from ", length(id), 
                 " annotated genomes ..", sep = ""), file = logfile)
      dna[id] <- sapply(gi[id], extractLocus, xml = xml, 
                        locus = x@locus)
    }
    
    ## checkpoint: was there failure to retrieve a single gene
    check <- which(dna == "")
    if ( length(check) > 0 ){
      stop("sequence retrieval failed for", paste("\n", sql.wrap(gi[check], term = "gi", BOOL = NULL)))
    }
      
    ## transform list in data frame
    ## ----------------------------
    if ( x@taxon@species.list ){
      taxon <- gsub(" ", "_", taxon)
    } else {
      taxon <- organism
    }
    seqs <- data.frame(gi = gi,
                       taxon = taxon,
                       spec_ncbi = organism,
                       status = "raw",
                       genom = topology,
                       npos = sapply(dna, nchar),
                       distbenchmark = NA,
                       dna = dna,
                       stringsAsFactors = FALSE)
    
    seqs <- unique(seqs) # yes, there are duplicated UIDs returned by eSearch!
    
    ## write into pgSQL database
    ## -------------------------
    if ( nrow(seqs) > 0 ) {
      conn <- dbconnect(x@db)
      present <- dbGetQuery(conn, paste("SELECT gi FROM", acc.tab))
      if ( nrow(present) > 0 ){
        id <- seqs$gi %in% present$gi
        slog("\n..", length(which(id)), "duplicates removed ..", file = logfile)
        seqs <- seqs[!id, ]
      }
      dbWriteTable(conn, acc.tab, seqs, row.names = FALSE, append = TRUE)
      slog("\n..", nrow(seqs), "sequences written to", acc.tab, "", file = logfile)  
      dbDisconnect(conn)
    }
    
    ## check if there are UIDs left to fetch
    ## -------------------------------------
    if ( nrow(seqs) == retmax){
      retstart <- i * retmax
      i <- i + 1
    } else {
      break
    } 
  } # END OF WHILE-loop
}