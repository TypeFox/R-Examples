dbReadReference <- function(dbc, gene){
  
  if ( !inherits(dbc, "PostgreSQLConnection") )
    stop("object 'dbc' is not a valid PostgreSQL connection")
  
  sql <- paste("SELECT rank, reference FROM reference WHERE",
               sql.wrap(gene, term = "gene"))
  sql <- dbGetQuery(dbc, sql)
  if ( nrow(sql) > 0 ){
    b <- as.DNAbin(strsplit(sql$reference, ""))
    names(b) <- sql$rank
  } else {
    b <- FALSE
  }
  b
}