
#' Vertex ids of a graph
#'
#' @param graph The graph.
#' @return Character vector of vertex ids.
#'
#' @export
#' @examples
#' G <- graph(list(A = c("B", "C"), B = "C", C = "A"))
#' vertex_ids(G)

vertex_ids <- function(graph) {
  graph <- as_graph_adjlist(graph)
  names(graph)
}

#' Vertices of a graph, with metadata
#'
#' @param graph The graph.
#' @return Character vector of vertex names.
#'
#' @family simple queries
#' @export
#' @examples
#' bridges <- graph(list(
#'   "Altstadt-Loebenicht" = c(
#'     "Kneiphof",
#'     "Kneiphof",
#'     "Lomse"
#'   ),
#'   "Kneiphof" = c(
#'     "Altstadt-Loebenicht",
#'     "Altstadt-Loebenicht",
#'     "Vorstadt-Haberberg",
#'     "Vorstadt-Haberberg",
#'     "Lomse"
#'   ),
#'   "Vorstadt-Haberberg" = c(
#'     "Kneiphof",
#'     "Kneiphof",
#'     "Lomse"
#'   ),
#'   "Lomse" = c(
#'     "Altstadt-Loebenicht",
#'     "Kneiphof",
#'     "Vorstadt-Haberberg"
#'   )
#' ))
#' vertices(bridges)

vertices <- function(graph) {
  graph <- as_graph_data_frame(graph)
  graph$nodes
}

#' Edges of a graph
#'
#' @param graph The graph
#' @return Data frame of edge data and metadata. The tail and head
#'   vertices are in the fist two columns. The rest of the columns are
#'   metadata.
#'
#' @family simple queries
#' @export
#' @examples
#' bridges <- graph(list(
#'   "Altstadt-Loebenicht" = c(
#'     "Kneiphof",
#'     "Kneiphof",
#'     "Lomse"
#'   ),
#'   "Kneiphof" = c(
#'     "Altstadt-Loebenicht",
#'     "Altstadt-Loebenicht",
#'     "Vorstadt-Haberberg",
#'     "Vorstadt-Haberberg",
#'     "Lomse"
#'   ),
#'   "Vorstadt-Haberberg" = c(
#'     "Kneiphof",
#'     "Kneiphof",
#'     "Lomse"
#'   ),
#'   "Lomse" = c(
#'     "Altstadt-Loebenicht",
#'     "Kneiphof",
#'     "Vorstadt-Haberberg"
#'   )
#' ))
#' edges(bridges)

edges <- function(graph) {
  graph <- as_graph_data_frame(graph)
  graph$edges
}

#' Order of a graph
#'
#' The order of the graph is the number of vertices.
#'
#' @param graph The graph.
#' @return Numeric scalar, the number of vertices.
#'
#' @family simple queries
#' @export
#' @examples
#' G <- graph(list(A = c("B", "C"), B = "C", C = "A"))
#' order(G)

order <- function(graph) {
  graph <- as_graph_adjlist(graph)
  length(graph)
}

#' The size of the graph is the number of edges
#'
#' @param graph The graph.
#' @return Numeric scalar, the number of edges.
#'
#' @export
#' @examples
#' G <- graph(list(A = c("B", "C"), B = "C", C = "A"))
#' size(G)

size <-function(graph) {
  graph <- as_graph_data_frame(graph)
  nrow(graph$edges)
}

#' Adjacent vertices for all vertices in a graph
#'
#' A vertex is adjacent is it is either a successor, or a predecessor.
#'
#' @param graph The graph.
#' @return A named list of character vectors, the adjacent vertices
#'   for each vertex.
#'
#' @family simple queries
#' @export
#' @examples
#' G <- graph(list(A = c("B", "C"), B = "C", C = "A"))
#' adjacent_vertices(G)

adjacent_vertices <- function(graph) {
   graph <- as_graph_adjlist(graph)
   merge_named_lists(graph, transpose(graph))
}

#' Predecessors and successors
#'
#' @param graph Input graph
#' @return Named list of character vectors, the predecessors or
#'   the successors of each vertex.
#'
#' @export
#' @examples
#' G <- graph(list(A = c("B", "C"), B = "C", C = "A"))
#' predecessors(G)
#' successors(G)

predecessors <- function(graph) {
  graph <- as_graph_adjlist(graph)
  unclass(transpose(graph))
}

#' @rdname predecessors
#' @export

successors <- function(graph) {
  unclass(as_graph_adjlist(graph))
}

#' Incident edges
#'
#' @param graph Input graph.
#' @param mode Whether to use \code{out} edges, \code{in} edges or
#'   \code{all} edges.
#' @return A list of data frames, each a set of edges.
#'
#' @export
#' @examples
#' G <- graph(list(A = c("B", "C"), B = "C", C = "A"))
#' incident_edges(G, mode = "out")
#' incident_edges(G, mode = "in")
#' incident_edges(G, mode = "all")

incident_edges <- function(graph, mode = c("out", "in", "all", "total")) {

  mode <- match.arg(mode)

  V <- vertex_ids(graph)
  edges <- edges(as_graph_data_frame(graph))

  res <- if (mode == "out") {
    lapply(
      V,
      function(w) drop_rownames(edges[w == edges[,1], ])
    )

  } else if (mode == "in") {
    lapply(
      V,
      function(w) drop_rownames(edges[w == edges[,2], ])
    )

  } else {
    lapply(
      V,
      function(w) drop_rownames(edges[w == edges[,1] | w == edges[,2], ])
    )
  }

  names(res) <- V
  res
}

#' Degree of vertices
#'
#' @param graph Input graph.
#' @param mode Whether to calculate \code{out}-degree, \code{in}-degree,
#'   or the \code{total} degree.
#' @return Named numeric vector of degrees.
#'
#' @export
#' @examples
#' G <- graph(list(A = c("B", "C"), B = "C", C = "A"))
#' degree(G, mode = "out")
#' degree(G, mode = "in")
#' degree(G, mode = "total")

degree <- function(graph, mode = c("out", "in", "total", "all")) {

  mode <- match.arg(mode)

  graph <- as_graph_adjlist(graph)
  V <- vertex_ids(graph)

  if (mode == "out") {
    vapply(graph, length, 1L)

  } else if (mode == "in") {
    vapply(transpose(graph), length, 1L)

  } else if (mode == "total" || mode == "all") {
    vapply(graph, length, 1L) + vapply(transpose(graph), length, 1L)
  }
}
