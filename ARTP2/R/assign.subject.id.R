
assign.subject.id <- function(null){
  
  if(!('SUBID' %in% colnames(null))){
    msg <- 'data frame \'data\' should have a column \'SUBID\''
    stop(msg)
  }
  
  id <- which(colnames(null) == 'SUBID')
  rownames(null) <- as.character(null$SUBID)
  null <- null[, -id, drop = FALSE]
  null
  
}
