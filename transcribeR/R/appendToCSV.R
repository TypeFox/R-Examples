appendToCSV <- function(existing.csv, row.frame, append = TRUE, sep=",", row.names=FALSE, col.names=FALSE) {
  if (is.null(existing.csv) || file.exists(existing.csv) == FALSE) {
    warning("The file was not found or path was wrong, \"file.name=" + existing.csv+ "\"")
    return(FALSE)
  }
  row.names(row.frame) <- NULL
  write.table(row.frame, 
              file = existing.csv,
              append=append,
              sep=sep, row.names=FALSE, 
              col.names = FALSE)
  return(TRUE)
}