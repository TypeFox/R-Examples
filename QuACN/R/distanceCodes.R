.distanceCodes <- function(g, dist=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(dist))
    dist <- distanceMatrix(g)

  maxdist <- max(dist)

  # for every row in the distance matrix,
  t(apply(dist, 1, function(row) {
    # count how often each value occurs
    result <- table(row, exclude=0)

    # and fill each row up with zeroes up to its diameter
    # in order to make apply() return a pretty matrix
    missing <- maxdist - length(result)
    if (missing > 0)
      result <- c(result,  rep(0, missing))

    names(result) <- 1:maxdist
    result
  }))
}
