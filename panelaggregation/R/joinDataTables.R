#' Joins two data.tables based on keys
#' 
#' This function joins two data.table objects, given a common key, which can have different names in the two tables.
#' In the latter case, the sequence of the names is crucial. Make sure that the key columns match exactly.
#' @author Matthias Bannert, Gabriel Bucur
#' @param dt_1 first data.table 
#' @param dt_2 second data.table
#' @param key_1 character vector of key columns for first data.table
#' @param key_2 character vector of key columns for second data.table
#' @return joined data.table
#' @export 
joinDataTables <- function(dt_1, dt_2, key_1, key_2 = key_1) {
  old_key1 <- key(dt_1)
  old_key2 <- key(dt_2)
  
  stopifnot(length(key_1) == length(key_2))
  
  setkeyv(dt_1, key_1)
  setkeyv(dt_2, key_2)
  
  res_dt <- dt_1[dt_2]
  
  setkeyv(dt_1, old_key1)
  setkeyv(dt_2, old_key2)
  
  res_dt
}
