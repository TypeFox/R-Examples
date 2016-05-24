#' Get basic network information as a string.
#'
#' @template arg_network
#' @param ... [any]\cr
#'   Not used at the moment.
#' @return [\code{character(1)}]
#' @export
as.character.Network = function(x, ...)   {
  n.points = getNumberOfNodes(x)
  n.clusters = getNumberOfClusters(x)

  char = if (!is.null(x$name)) paste0(x$name, "\n") else ""
  char = paste0(char, "#Nodes: ", n.points)
  if (n.clusters > 1L)
    char = paste0(char, ", #Clusters: ", n.clusters)
  if (hasAttributes(x, "morphed"))
    char = paste0(char, "(Morphing coefficient ", attr(x, "morphing.grade"), ")")
  char
}
