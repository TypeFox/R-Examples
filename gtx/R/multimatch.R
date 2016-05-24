multimatch <- function(query, target, values, sep = ",", use.unique = TRUE) {
  stopifnot(length(values) == length(target))
  if (use.unique) return(sapply(1:length(query), function(ii) return(paste(unique(values[target == query[ii]]), collapse = sep))))
  if (!use.unique) return(sapply(1:length(query), function(ii) return(paste(values[target == query[ii]], collapse = sep))))
  return(invisible(NULL)) # this should never happen
}  

  
