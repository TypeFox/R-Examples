
#' Is the graph weighted?
#'
#' @param graph The graph.
#'
#' @export
#' @examples
#' G <- graph(
#'   data.frame(
#'     stringsAsFactors = FALSE,
#'     id = c("a", "b", "c", "d")
#'   ),
#'   data.frame(
#'     stringsAsFactors = FALSE,
#'     from   = c("a", "a", "b", "b", "c"),
#'     to     = c("b", "d", "d", "c", "a"),
#'     weight = c( 1 ,  2 ,  1 ,  3 ,  2 )
#'   )
#' )
#' is_weighted(G)
#'
#' G2 <- graph(
#'   data.frame(
#'     stringsAsFactors = FALSE,
#'     id = c("a", "b", "c", "d")
#'   ),
#'   data.frame(
#'     stringsAsFactors = FALSE,
#'     from   = c("a", "a", "b", "b", "c"),
#'     to     = c("b", "d", "d", "c", "a")
#'   )
#' )
#' is_weighted(G2)

is_weighted <- function(graph)
  UseMethod("is_weighted")

#' @method is_weighted simplegraph_adjlist
#' @export

is_weighted.simplegraph_adjlist <- function(graph) {
  FALSE
}

#' @method is_weighted simplegraph_df
#' @export

is_weighted.simplegraph_df <- function(graph) {
  "weight" %in% colnames(graph$edges)
}

#' Vertex strength: sum of weights of incident edges
#'
#' This is also called weighed degree.
#'
#' For non-weighted graphs, the degree is returned as a
#' fallback.
#'
#' @param graph Input graph.
#' @param mode Whether to consider incoming (\code{in}),
#'   outgoing (\code{out}) or all (\code{total}) edges.
#' @return Named numeric vector.
#'
#' @export
#' @examples
#' G <- graph(
#'   data.frame(
#'     stringsAsFactors = FALSE,
#'     id = c("a", "b", "c", "d")
#'   ),
#'   data.frame(
#'     stringsAsFactors = FALSE,
#'     from   = c("a", "a", "b", "b", "c"),
#'     to     = c("b", "d", "d", "c", "a"),
#'     weight = c( 1 ,  2 ,  1 ,  3 ,  2 )
#'   )
#' )
#' strength(G)
#'
#' G2 <- graph(
#'   data.frame(
#'     stringsAsFactors = FALSE,
#'     id = c("a", "b", "c", "d")
#'   ),
#'   data.frame(
#'     stringsAsFactors = FALSE,
#'     from   = c("a", "a", "b", "b", "c"),
#'     to     = c("b", "d", "d", "c", "a")
#'   )
#' )
#' strength(G2)

strength <- function(graph, mode = c("out", "in", "total", "all")) {

  if (!is_weighted(graph)) return(degree(graph = graph, mode = mode))

  mode <- match.arg(mode)

  graph <- as_graph_data_frame(graph)

  inc <- incident_edges(graph, mode = mode)
  vapply(inc, function(x) sum(x$weight), 1.0)
}
