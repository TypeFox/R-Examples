#' @export
print.Network = function(x, ...) {
  clustered = testClass(x, "ClusteredNetwork")
  type = if (clustered) "Clustered" else "Simple"
  catf("%s %i-dimensional network.", type, ncol(x$coordinates))
  catf("Name:               %s", x$name)
  catf("Comment(s):         %s", collapse(x$comment, "\n"))
  catf("Edge weight type:   %s", x$edge.weight.type)
  catf("Number of nodes:    %i", getNumberOfNodes(x))
  if (clustered) {
    catf("Number of clusters: %i", getNumberOfClusters(x))
  }
  if (hasDepots(x)) {
    catf("Number of depots:   %i", getNumberOfDepots(x))
  }
  catf("Head of coordinate matrix:")
  x$coordinates = round(x$coordinates, digits = 2L)
  df = as.data.frame(x, include.extra = TRUE)
  print(head(df, n = 5L))
  catf("...")
}
