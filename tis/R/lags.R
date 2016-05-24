lags <- function(x, lags, name = ""){
  ## returns a multivariate tis series of specified lags of x with appropriate
  ## colnames 
  z <- 0
  for(i in lags) z <- cbind(z, lag(x, i))
  z <- z[,-1, drop = F]

  lnames <- paste("(+", lags, ")", sep = "")
  lagNames <- gsub("\\+0", "0", gsub("\\+-", "-", lnames))
                       
  if(!missing(name) && length(name) == NCOL(x)) cn <- name
  else if(is.null(cn <- colnames(x))) cn <- letters[1:NCOL(x)]
  
  colnames(z) <- as.vector(outer(cn, lagNames, paste, sep = ""))
  z
}

Lags <- function(x, lags, name = "") lags(x, -lags, name)
