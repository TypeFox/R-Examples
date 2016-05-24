dbWriteDNA <- function(conn, tab.name, seqs, enforce.binomial = TRUE, status){
  
  if ( is.matrix(seqs) ) seqs <- as.list(seqs)
  seqs <- as.character(seqs)
  len <- sapply(seqs, length)
  seqs <- lapply(seqs, c2s)
  
  id <- strsplit(names(seqs), "_")
  id <- do.call(rbind, id)
  ## will cut subgeneric names:
  gi = id[, ncol(id)] 
  taxon = id[, -ncol(id), drop = FALSE]
  if ( enforce.binomial) taxon <- taxon[, 1:2]
  taxon <- apply(taxon, 1, paste, collapse = "_")
  
  status <- ifelse ( missing(status), "",  paste("status='", status, "', ", sep = ""))
    
  for ( i in seq_along(seqs) ){
    SQL <- paste("UPDATE ", tab.name, " SET ", status, "npos='", len[i], "', dna='", 
                 seqs[[i]], "'", " WHERE gi='", gi[i], "' AND taxon='", 
                 taxon[i], "';", sep = "")
    dbSendQuery(conn, SQL)
  }
}
