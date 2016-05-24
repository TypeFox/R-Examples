is.what <- function(object, verbose=FALSE)
{
  do.test <- function(test, object)
  {
    result <- try(get(test)(object), silent=TRUE)
    if(!is.logical(result) || is.na(result) || length(result)!=1)
      result <- NULL
    return(result)
  }

  ## Get all names starting with "is."
  is.names <- unlist(sapply(search(), function(name) ls(name,pattern="^is\\.")))

  ## Narrow to functions
  is.functions <- is.names[sapply(is.names, function(x) is.function(get(x)))]

  tests <- sort(unique(is.functions[is.functions!="is.what"]))
  results <- suppressWarnings(unlist(sapply(tests, do.test, object=object)))

  if(verbose)
    output <- data.frame(is=ifelse(results,"T","."))
  else
    output <- names(results)[results]

  return(output)
}
