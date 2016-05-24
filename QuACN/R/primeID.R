primeID <- function(g, deg=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(deg))
    deg <- graph::degree(g)

  n <- nodes(g)
  e <- edges(g)
  .weightedPathSum(e, NULL, function(i, from, to) {
    pair <- sort(deg[c(from, to)])
    prime <- if (identical(pair, c(1, 2)))        2
             else if (identical(pair, c(1, 3)))   3
             else if (identical(pair, c(1, 4)))   5
             else if (identical(pair, c(2, 2)))   7
             else if (identical(pair, c(2, 3)))   11
             else if (identical(pair, c(2, 4)))   13
             else if (identical(pair, c(3, 3)))   17
             else if (identical(pair, c(3, 4)))   19
             else                                 23
    1 / sqrt(prime)
  })
}
