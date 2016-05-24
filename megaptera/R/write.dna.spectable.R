write.dna.spectable <- function(conn, tab.name, dna, overwrite = FALSE){
  
  ## prepare data
  ## ------------
  if ( is.matrix(dna) ) dna <- as.list(dna)
  dna <- as.character(dna)
  npos <- sapply(dna, length)
  dna <- lapply(dna, c2s)
  spec <- names(dna)
  
  spec.in.db <- dbGetQuery(conn, paste("SELECT spec FROM", tab.name))$spec
  insert <- rep(TRUE, length(spec))
  insert[spec %in% spec.in.db] <- FALSE
  
  for ( i in seq_along(dna) ){
    SQL <- ifelse(insert[i],
                  paste("INSERT INTO", tab.name, 
                        "(spec, block, npos, dna)",  
                        "VALUES (", 
                        sql.wrap(spec[i], term = NULL), ",",
                        "'raw',",
                        sql.wrap(npos[i], term = NULL), ",", 
                        sql.wrap(dna[[i]], term = NULL), ")"),
                  paste("UPDATE", tab.name, 
                        "SET block = 'raw',", 
                        sql.wrap(npos[i], term = "npos"), ",", 
                        sql.wrap(dna[[i]], term = "dna"), 
                        "WHERE", sql.wrap(spec[i])))
    dbSendQuery(conn, SQL)
  }
}
