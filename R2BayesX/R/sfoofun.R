sfoofun <- function(x, xt = NULL, ...)
{
  if(is.null(x) || is.null(xt))
    return(NULL)
  x <- rval <- deparse(substitute(xt), backtick = TRUE, width.cutoff = 500L)
  if(is.list(xt) && length(xt)>1)
    for(i in 1L:length(xt))
      if(inherits(xt[[i]], "bnd") || inherits(xt[[i]], "gra") || inherits(xt[[i]], "list")) {
        rval <- strsplit(x, ",", " ")[[1L]]
        if(length(rval) > 1L)
          rval <- rval[i]
      }
  rval <- splitme(rval)
  if(length(rval) > 5L)
    if(resplit(rval[1L:5L]) == "list(")
      rval <- rval[6L:length(rval)]
  if(rval[length(rval)] == ")")
    rval <- rval[1L:(length(rval) - 1)]
  if(any(grepl("=", rval)))
    rval <- rval[(grep("=", rval) + 2L):length(rval)]
  rval <- resplit(rval)
   
  return(rval)
}

