rrd <-
function(x) 
{
  xo <- x
  x <- splitme(x)
  go <- TRUE
  rval <- NULL
  for(i in 1L:length(x)) {
    if(x[i] == ",")
      go <- FALSE
    if(go)
      rval <- c(rval, x[i])
  }
  rval <- resplit(c(rval, ")"))
  rval <- gsub("))", ")", rval)
  if(grepl("):", xo, fixed = TRUE))
    rval <- paste(rval, strsplit(xo, "):", fixed = TRUE)[[1L]][2L], sep = ":")
  
  return(rval)
}

