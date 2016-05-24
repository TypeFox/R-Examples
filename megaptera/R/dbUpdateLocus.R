# SEARCH AND DOWNLOAD SEQUENCES in parallel execution
# package: megaptera
# called by: 
# author: Christoph Heibl
# last change 2014-08-20

dbUpdateLocus <- function(x, maintain = FALSE){
  
  ## DEFINITIONS
  ## -----------
  gene <- x@locus@sql
  acc.tab <- paste("acc", gsub("^_", "", gene), sep = "_")
  cols <- paste(gene, c("gb", "sel", "blocks"), sep = "_")
  
  ## open database connection
  ## ------------------------
  conn <- dbconnect(x@db)
 
  ## read taxonomy table
  ## -------------------
  tax <- dbReadTaxonomy(conn)[, "spec", drop = FALSE]
  
  ## set up (if it does not exist) and read locus table
  ## --------------------------------------------------
  if ( !dbExistsTable(conn, "locus") ){
    dbWriteTable(conn, "locus", tax, row.names = FALSE)
    dbSendQuery(conn, "ALTER TABLE locus ADD CONSTRAINT locus_pk PRIMARY KEY (spec)")
  }
  locus <- dbReadTable(conn, "locus")
  
  ## add new rows = species (if missing)
  ## -----------------------------------
  specMissing <- setdiff(tax$spec, locus$spec)
  if ( length(specMissing) > 0 ) {
    SQL <- paste("INSERT INTO locus (spec) VALUES (", 
                 sql.wrap(specMissing, term = NULL, BOOL = NULL), ")", 
                 sep = "")
    lapply(SQL, dbSendQuery, conn = conn)
  }
  
  ## add new columns = loci (if missing)
  ## -----------------------------------
  if ( length(grep(cols[1], names(locus))) == 0 ) {
    SQL <- paste("ALTER TABLE locus ADD COLUMN", 
                 paste(cols, c(rep("integer", 2), "text")))
    lapply(SQL, dbSendQuery, conn = conn)
  }
  
  if ( maintain ){
    dbDisconnect(conn)
    return(NULL)
  } 
  ## update locus table: how many sequences per species on NCBI?
  ## -----------------------------------------------------------
  tax.found <- paste("SELECT taxon FROM", acc.tab)
  tax.found <- dbGetQuery(conn, tax.found)
  tax.found <- table(tax.found$taxon) ## count number of accesssion
  
  SQL <- paste("UPDATE locus SET",  sql.wrap(0, term = cols[1]))
  if ( nrow(tax.found) > 0 ){
    SQL <- c(SQL, 
             paste("UPDATE locus",
                   "SET", sql.wrap(tax.found, term = cols[1], BOOL = NULL),
                   "WHERE", sql.wrap(names(tax.found), BOOL = NULL)))
  }
  lapply(SQL, dbSendQuery, conn = conn)
  
  ## close database connection
  ## ------------------------
  dbDisconnect(conn)
}