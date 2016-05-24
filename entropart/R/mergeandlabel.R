mergeandlabel <-
function (x, y) {
  # Merge two dataframes by row names, keep all rows
  Merged <- merge(x, y, by=0, all=TRUE)
  # Row names are in the first column. Name Rows and eliminate the column
  row.names(Merged) <- Merged[, 1]
  Merged <- Merged[, -1]
  # Replace NA's by zeros
  Merged[is.na(Merged)] <- 0
  # Set unique names to avoid warnings when Reduce() is used
  names(Merged) <- paste("C", 1:ncol(Merged), sep="")
  return(Merged)
}
