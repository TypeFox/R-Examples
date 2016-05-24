
#' @method graph data.frame
#' @export

graph.data.frame <- function(x, y, ...) {
  sanitize(
    structure(
      list(nodes = x, edges = y),
      class = c("simplegraph_df", "simplegraph")
    )
  )
}

#' @method sanitize simplegraph_df
#' @export

sanitize.simplegraph_df <- function(x, ...) {
  if (!"nodes" %in% names(x)) stop("No vertices found in data frame graph")
  if (ncol(x$nodes) == 0) stop("Vertex data frame must contain ids")
  if (!is.character(x$nodes[,1])) {
    stop("First column must contain (character) ids in vertex data frame")
  }
  if (!"edges" %in% names(x)) stop("No edges found in data frame graph")
  if (ncol(x$edges) < 2 ||
      !is.character(x$edges[,1]) ||
      !is.character(x$edges[,2])) {
    stop("First two columns must contain (character) ids in edge data frame")
  }
  if (any(! x$edges[,1] %in% x$nodes[,1]) ||
      any(! x$edges[,2] %in% x$nodes[,1])) {
    stop("Unknown vertex id in edge data frame")
  }
  x
}

as_graph_data_frame <- function(x, ...)
  UseMethod("as_graph_data_frame")

as_graph_data_frame.simplegraph_df <- function(x, ...) {
  x
}

as_graph_data_frame.simplegraph_adjlist <- function(x, ...) {
  nodes <- data.frame(
    stringsAsFactors = FALSE,
    name = names(x)
  )
  edges <- data.frame(
    stringsAsFactors = FALSE,
    row.names = NULL,
    from = rep(names(x), vapply(x, length, 1L)),
    to = unlist(x)
  )
  graph(nodes, edges)
}
