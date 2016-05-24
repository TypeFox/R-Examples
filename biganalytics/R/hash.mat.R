hash.mat <- function(x, col)
{
  col <- bigmemory:::cleanupcols(col, ncol(x), colnames(x))
  if (colmin(x, col)<1)
    stop("Error: minimum value in specified column should be 1 or more.")
  return(matrix(.Call('MatrixHashRanges', x@address, as.double(col), 
                      PACKAGE="biganalytics"), ncol=2, byrow=TRUE))
}


