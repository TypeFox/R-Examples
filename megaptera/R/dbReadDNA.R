# READ SEQUENCES from DATABASE
# package: megaptera
# author: Christoph Heibl
# last update: 2014-09-05

dbReadDNA <- function(conn, tab.name, taxon, regex = FALSE, max.bp, max.dist,
                      enforce.binomial = FALSE, ignore.excluded = TRUE){
  
  if ( missing(taxon) ) taxon <- ".+"
  otaxon <- taxon
  
  ## esacape metacharacters in taxon names
  ## -------------------------------------
  if ( !regex ){
    taxon <- gsub(" ", "_", taxon)
    taxon <- gsub("[.]$", "[.]", taxon)
    taxon <- gsub("([(]|[+]|[)])", "[\\1]", taxon)
    taxon <- gsub("'", ".", taxon) # e.g.Acorus_sp._'Mt._Emei'
  }
  
  ## field names in <tab.name>
  ## -------------------------
  cols <- paste("SELECT column_name FROM information_schema.columns WHERE", 
                sql.wrap(tab.name, term = "table_name"), 
                "ORDER by ordinal_position")
  cols <- dbGetQuery(conn, cols)$column_name
  
  if ( "taxon" %in% cols ){
    ## retrieve sequences from acc.table
    ## ---------------------------------
    if ( enforce.binomial ) taxon <- paste(taxon, "$", sep = "")
    SQL <- paste("SELECT taxon, gi, dna FROM ", tab.name, 
                 " WHERE taxon ~ '^", taxon, "'", sep = "")
    if ( !missing(max.bp) ) SQL <- paste(SQL, "AND npos <=", max.bp)
    if ( !missing(max.dist) ) SQL <- paste(SQL, "AND distreference <=", max.dist)
    if ( ignore.excluded ) SQL <- paste(SQL, "AND status !~ 'excluded'")
#     print(SQL)
    seqs <- dbGetQuery(conn, SQL)
    if ( nrow(seqs) == 0) stop("no sequences for '", otaxon, "'")
    seq.names <- paste(seqs$taxon, seqs$gi, sep = "_")
  } else {
    ## retrieve sequences from spec.table
    ## ----------------------------------
    SQL <- ifelse(ignore.excluded,
                   " AND block !~ 'excluded'", "") 
    SQL <- paste("SELECT spec, dna FROM ", tab.name, 
                 " WHERE spec ~ '^", taxon, "'", SQL, sep = "")
    if ( !missing(max.bp) ) SQL <- paste(SQL, "AND npos <=", max.bp)
    seqs <- dbGetQuery(conn, SQL)
    if ( nrow(seqs) == 0 ) stop("no sequences for '", otaxon, "'")
    seq.names <- seqs$spec
    
  }
  seqs <- as.list(seqs$dna)
  seqs <- lapply(seqs, s2c)
  names(seqs) <- seq.names
  seqs <- as.DNAbin(seqs)
  
  ## convert to matrix if possible
  ## -----------------------------
  if ( length(unique(sapply(seqs, length))) == 1 ){
    seqs <- as.matrix(seqs)
    seqs <- deleteEmptyCells(seqs, quiet = TRUE)
  }
  
  return(seqs)  
}
