eventsJson2RMatrix = function(
  jsonData, 
  select=TRUE  # If a vector of column names, only specified properties are selected. 
) {
  data = jsonlite::fromJSON(paste("[", paste(jsonData, collapse=","), "]"))
  
  if (!is.logical(select)) {
    notFound = !select %in% colnames(data$properties)
    if (any(notFound)) {
      paste("Warning: some properties not found:", paste(select[notFound], collapse=", "))
      select = select[!notFound]
    }
  }
  
  r1 = data$event
  r2 = data$properties[, select, drop=F]
  
  d = matrix(NA, length(r1), 1 + ncol(r2))
  colnames(d) = c("event", colnames(r2))
  d[, 1] = as.character(r1)
  
  ## Fill empty array. Replace NULLs by NAs. Collapse array user properties.
  for (i in 1:ncol(r2)) {
    ri = r2[, i]
    
    ## Dictionaries look like arrays...
    if (!is.null(dim(ri)) && length(dim(ri) == 2)) {
      if (ncol(ri) == 0)
        d[, 1+i] = NA
      else
        d[, 1+i] = apply(ri, 1, function(z) paste(names(z), "=", z, sep="", collapse=",")) # Collapse dict.
    } else
      d[, 1+i] = unlist(lapply(ri, function(z) if(is.null(z)) NA else paste(z, collapse=","))) # Collapse array
  }
  d
}