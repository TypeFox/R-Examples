df2m <- function(x)
{
  if(!is.null(x)) {
    xattr <- attributes(x)
    nxa <- names(xattr)
    x$intnr <- x$paramnr <- x$varname <- NULL
    cn <- colnames(x)
    x <- as.matrix(x)
    rownames(x) <- 1L:nrow(x)
    colnames(x) <- rep(cn, length.out = ncol(x))
    for(k in 1L:length(nxa)) 
      if(all(nxa[k] != c("dim", "dimnames", "class", "names", "row.names"))) {
        attr(x, nxa[k]) <- xattr[[k]]
      }
    }

  return(x)
}

