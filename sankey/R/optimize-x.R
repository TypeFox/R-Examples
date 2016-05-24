
#' @importFrom simplegraph graph topological_sort order vertex_ids

optimize_x <- function(nodes, edges) {

  ## `simplegraph` object
  sgraph <- graph(nodes, edges)

  ## Reverse adjacency list of the graph
  adj <- predecessors(sgraph)

  levels <- structure(rep(-1, order(sgraph)), names = vertex_ids(sgraph))

  order <- topological_sort(sgraph)

  for (n in order) {
    pred_levels <- levels[ adj[[n]] ]
    levels[[n]] <- if (length(pred_levels) == 0) 0 else max(pred_levels) + 1
  }

  levels
}
