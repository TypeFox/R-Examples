dbDeleteSpec <- function(conn, spec, gene){
  
  if ( class(conn) == "sproj" ){
    if ( missing(gene) ) {
      gene <- conn$gene$sql
      all.genes <- FALSE
    } else {
      all.genes <- TRUE
    }
    conn <- dbConnect(PostgreSQL(), host = x$db[["host"]], 
                      port = x$db[["port"]], dbname = x$db[["dbname"]], 
                      user = x$db[["user"]], password = x$db[["password"]])
    
  } else {
    if ( missing(gene) ) all.genes <- TRUE
  }
  
  ## species
  spec <- gsub(" ", "_", spec)
  SQL <- sql.wrap(spec)
  
  if ( all.genes ){
    
    x <- dbGetQuery(conn, paste("SELECT * FROM taxonomy WHERE", SQL))
    cols <- grep("_gb", colnames(x))[-1]
    cols <- colnames(x)[cols][x[, cols] > 0]
    cols <- gsub("_gb", "", cols)
    SQL <- c(paste("DELETE FROM taxonomy WHERE", SQL),
             paste("DELETE FROM", cols, "WHERE", sql.wrap(spec, term = "taxon")))
    lapply(SQL, dbSendQuery, conn = conn)
  } else {
    ## all.genes = FALSE
    SQL <- sql.wrap(spec)
    cols <- paste(gene, c("gb", "sel", "blocks"), sep = "_")
    cols <- paste(cols, "NULL", sep = "=")
    cols <- paste(cols, collapse = ", ")
    dbSendQuery(conn, paste("UPDATE taxonomy SET", cols, "WHERE", SQL))
    dbSendQuery(conn, paste("DELETE FROM", gene, "WHERE", sql.wrap(spec, term = "taxon")))
  }
}