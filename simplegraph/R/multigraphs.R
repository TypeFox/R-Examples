
#' Is this a simple graph?
#'
#' A simple graph contains no loop and multiple edges.
#'
#' @param graph The input graph.
#' @return Logical scalar.
#'
#' @family multigraphs
#' @export
#' @examples
#' G <- graph(list(A = c("A", "B", "B"), B = c("A", "C"), C = "A"))
#' is_simple(G)
#'
#' G2 <- simplify(G)
#' is_simple(G2)

is_simple <- function(graph) {
  ! is_loopy(graph) && ! is_multigraph(graph)
}

#' Is this a loopy graph?
#'
#' A loopy graph has at least one loop edge: an edge from a vertex to
#' itself.
#'
#' @param graph The input graph.
#' @return Logical scalar.
#'
#' @family multigraphs
#' @export
#' @examples
#' G <- graph(list(A = c("A", "B", "B"), B = c("A", "C"), C = "A"))
#' is_loopy(G)
#'
#' G2 <- simplify(G)
#' is_loopy(G2)

is_loopy <- function(graph) {

  graph <- as_graph_adjlist(graph)

  for (n in names(graph)) {
    if (n %in% graph[[n]]) return(TRUE)
  }

  FALSE
}


#' Is this a multigraph?
#'
#' A multigraph has at least one pair or multiple edges,
#' edges connecting the same (ordered) pair of vertices.
#'
#' @param graph Input graph.
#' @return Logical scalar.
#'
#' @family multigraphs
#' @export
#' @examples
#' G <- graph(list(A = c("A", "B", "B"), B = c("A", "C"), C = "A"))
#' is_multigraph(G)
#'
#' G2 <- simplify(G)
#' is_multigraph(G2)

is_multigraph <- function(graph) {

  graph <- as_graph_adjlist(graph)

  for (adj in graph) {
    if (any(duplicated(adj))) return(TRUE)
  }

  FALSE
}

#' Remove multiple and loop edges from a graph
#'
#' @param graph Input graph.
#' @return Another graph, with the multiple and loop edges removed.
#'
#' @family multigraphs
#' @export
#' @examples
#' G <- graph(list(A = c("A", "B", "B"), B = c("A", "C"), C = "A"))
#' is_simple(G)
#'
#' G2 <- simplify(G)
#' is_simple(G2)

simplify <- function(graph) {
  remove_loops(remove_multiple(graph))
}

#' Remove loop edges from a graph
#'
#' @param graph Input graph
#' @return Graph, with loop edges removed.
#'
#' @family multigraphs
#' @export
#' @examples
#' G <- graph(list(A = c("A", "B", "B"), B = c("A", "C"), C = "A"))
#' is_loopy(G)
#' is_loopy(remove_loops(G))

remove_loops <- function(graph)
  UseMethod("remove_loops")

#' @method remove_loops simplegraph_df
#' @export

remove_loops.simplegraph_df <- function(graph) {
  graph$edges <- graph$edges[graph$edges[,1] != graph$edges[,2], ]
  row.names(graph$edges) <- seq_len(nrow(graph$edges))
  graph
}

#' @method remove_loops simplegraph_adjlist
#' @export

remove_loops.simplegraph_adjlist <- function(graph) {
  graph(
    structure(
      lapply(
        names(graph),
        function(n) setdiff(graph[[n]], n)
      ),
      names = names(graph)
    )
  )
}

#' Remove multiple edges from a graph
#'
#' @param graph Input graph.
#' @return Graph, without the multiple edges. (More precisely, from
#'   each set of multiple edges, only one, the first one, is kept.)
#'
#' @family multigraphs
#' @export
#' @examples
#' G <- graph(list(A = c("A", "B", "B"), B = c("A", "C"), C = "A"))
#' is_multigraph(G)
#' is_multigraph(remove_multiple(G))

remove_multiple <- function(graph)
  UseMethod("remove_multiple")

#' @method remove_multiple simplegraph_df
#' @export

remove_multiple.simplegraph_df <- function(graph) {
  graph$edges <- graph$edges[ ! duplicated(graph$edges[,1:2]), ]
  graph
}

#' @method remove_multiple simplegraph_adjlist
#' @export

remove_multiple.simplegraph_adjlist <- function(graph) {
  graph(lapply(graph, unique))
}
