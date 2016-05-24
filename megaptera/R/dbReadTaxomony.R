dbReadTaxonomy <- function(conn, subset){
  
  if ( !inherits(conn, "PostgreSQLConnection") )
    stop("object 'conn' is not a valid PostgreSQL connection")
  
  if ( !dbExistsTable(conn, "taxonomy") )
    stop("no taxonomy table - see ?dbUpdateTaxonomy for help")
  
  tax <- dbReadTable(conn, "taxonomy")
  
  ## subsetting taxonomy ..
  ## ----------------------
  if ( !missing(subset) ){
    ## .. based on sequence names or ..
    ## --------------------------------
    if ( inherits(subset, "DNAbin") ){
      if ( is.list(subset) ) sset <- names(subset)
      if ( is.matrix(subset) ) sset <- rownames(subset)
      tax <- tax[tax$spec %in% sset, ]
    }
    ## .. based on sequence names or ..
    ## --------------------------------
    if ( inherits(subset, "phylo") ){
      tax <- tax[tax$spec %in% subset$tip.label, ]
    }
    ## .. based on <spec.*>
    ## --------------------
    if ( is.character(subset) ){
      subset <- paste("SELECT spec FROM", subset, 
                      "WHERE block !~ 'excluded'")
      subset <- dbGetQuery(conn, subset)$spec
      tax <- tax[tax$spec %in% subset, ]
    }
  }
  
  ## remove trailing white spaces
  ## happens often when people prepare taxon lists in Excel
  ## ------------------------------------------------------
  tws <- grep(" $|_$", tax[, "spec"])
  if ( length(tws) > 0 ){
    tax[, "spec"] <- gsub(" $|_$", "", tax[, "spec"])
    warning("trailing white space removed in", 
            paste("\n - ", head(tax[, "spec"][tws], 3), sep = ""), 
            "\n - and ", length(tws) - 3, " more species")    
  }
  
  tax <- sqlTaxonomyHeader(tax) # should be obsolete
  tax
}