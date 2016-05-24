mixNamedVec <-
function(x, y, resolve.conflicts = c("x", "y"), warn = TRUE) {
  if (!(is.vector(x) & is.vector(y))) stop("x and y must be vectors")
  resolve.conflicts <- match.arg(resolve.conflicts)
  
  xnames <- names(x)
  ynames <- names(y)
  
  duplicates <- xnames[xnames %in% ynames]
  
  if ((warn) & !is.null(duplicates)) 
    warning(paste("name duplicates found:", paste(dQuote(duplicates), collapse=", ")))
  if (resolve.conflicts == "x")
    return(c(y[!(ynames %in% xnames)], x))
  if (resolve.conflicts == "y")
    return(c(x[!(xnames %in% ynames)], y))
}

