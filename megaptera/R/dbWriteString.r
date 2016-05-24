dbWriteString <- function(conn, name, value){
  
  if ( is.matrix(value) ) value <- as.list(value)
  
  for ( i in seq_along(value) ){
    SQL <- paste("INSERT INTO ", name, 
                 " (spec, dna) VALUES ('", names(value)[i],
                 "', '", 
                 paste(as.character(value[[i]]), collapse = ""),
                 "')", 
                 sep = "")
    dbSendQuery(conn, SQL)
  }
}