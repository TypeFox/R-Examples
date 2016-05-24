# updated 2014-01-30
check.Species <- function(conn, spec, gene){
  
  close.after <- FALSE
  if ( class(conn) == "sproj" ){
    gene <- conn$gene
    conn <- dbConnect(PostgreSQL(), host = x$db[["host"]], 
                      port = x$db[["port"]], dbname = x$db[["dbname"]], 
                      user = x$db[["user"]], password = x$db[["password"]])
    close.after <- TRUE
  }
  
  ## taxonomy
  ## --------
  sp <- gsub(" ", "_", spec)
  mk <- sql.conform(gene)
  tab <- dbGetQuery(conn, paste("SELECT * FROM taxonomy WHERE", 
                                sql.wrap(sp)))
  mk <- tab[, grep(mk, names(tab))]
  tab <- tab[, -grep("_gb|_sel|_blocks", names(tab))]
  
  ## assemble output object
  x <- list()
  x$syn <- tab$synonmys
  tab <- tab[, -grep("synonyms", names(tab))]
  x$tax <- tab[1, , drop = TRUE]
  x$mk <- mk
  x$genbank <-dbGetQuery(conn, paste("SELECT gi, genom, distbenchmark", 
                                     "FROM", sql.conform(gene), 
                                     "WHERE", sql.wrap(sp, term = "taxon")))
  if ( close.after ) dbDisconnect(conn)
  ## output
  cat("--- Taxonomy ---\n")
  cat(paste(names(x$tax), ": ", x$tax, "\n", sep = ""))
  cat("\n--- GenBank ---\n")
  nseq <- nrow(x$genbank)
  if ( nseq == 0 ){
    cat("No sequences found at GenBank")
  } else {
    cat(paste("Number of sequences: ", nseq, "\n", sep = ""))
    cat(paste("Distance from reference sequence: ", paste(round(mean(x$genbank$distbenchmark, na.rm = TRUE), 5), 
                                                          " (", 
                                                          paste(range(x$genbank$distbenchmark, na.rm = TRUE), 
                                                                collapse = "-"), 
                                                          ")", sep = ""), 
              "\n", sep = ""))
    cat("\n")
    cat(paste("Sequences selected for contig: ", x$mk[1, 2], "\n", sep = ""))
    cat("\n")
    cat(paste("Sequence to find in: ", x$mk[1, 3], "\n", sep = ""))
  }
}