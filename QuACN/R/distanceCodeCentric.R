distanceCodeCentric <- function(g, dist=NULL) {

  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dist))
    dist <- distanceMatrix(g)

  n <- numNodes(g)

  # calculate frequency of rows in distance matrix
  # by first pasting the columns into strings
  dist <- apply(dist, 1, function(row) do.call(paste, as.list(sort(row))))
  row_freq <- table(dist)
  row_freq <- row_freq[row_freq != 0]
  p <- row_freq / n

  -sum(p * log2(p))
}
