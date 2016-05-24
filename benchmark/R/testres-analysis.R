

contsigseq <- function(x, significance) {
  #x <- subset(x, subset = tests == "pvalue", select = value)
  x <- x[x$tests == "pvalue", "value"]
  x <- x < significance

  i <- length(x)
  while ( x[i] == TRUE )
    i <- i - 1

  if ( i == length(x) )
    return(NA)

  i + 1
}


