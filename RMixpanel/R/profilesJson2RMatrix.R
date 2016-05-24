profilesJson2RMatrix <- function(
  jsonData, 
  select=TRUE  # If a vector of column names, only specified columns are selected. 
) {
  if (!is.logical(select)) {
    notFound = !select %in% colnames(jsonData$results[, 2])
    if (any(notFound)) {
      paste("Warning: some properties not found:", paste(select[notFound], collapse=", "))
      select = select[!notFound]
    }
  }
  
  r1 = jsonData$results[, 1]
  r2 = jsonData$results[, 2][, select, drop=F]
  
  d = matrix(NA, length(r1), 1 + ncol(r2))
  colnames(d) = c("distinct_id", colnames(r2))
  d[, 1] = as.character(r1)
  
  ## Fill empty array. Replace NULLs by NAs. Collapse array profile properties.
  for (i in 1:ncol(r2))
    d[, 1+i] = unlist(lapply(r2[, i], function(z) if(is.null(z)) NA else paste(z, collapse=",")))
  d
}
