
#' @method graph list
#' @export

graph.list <- function(x, ...) {
  sanitize(
    structure(
      x,
      class = c("simplegraph_adjlist", "simplegraph", class(x))
    )
  )
}

#' @method sanitize simplegraph_adjlist
#' @export

sanitize.simplegraph_adjlist <- function(x, ...) {
  if (is.null(names(x))) stop("Adjacency list must be named")
  if (any(names(x) == "")) stop("Names must be non-empty in adjacency list")
  if (any(duplicated(names(x)))) stop("Duplicated names in adjacency list")
  if (any(!vapply(x, is.character, TRUE))) {
    stop("Adjacency list must contain character vectors")
  }
  if (any(!vapply(x, function(l) all(l %in% names(x)), TRUE))) {
    stop("Unknown vertices in adjacency list")
  }
  x
}

as_graph_adjlist <- function(x, ...)
  UseMethod("as_graph_adjlist")

as_graph_adjlist.simplegraph_adjlist <- function(x, ...) {
  x
}

#' @importFrom utils modifyList

as_graph_adjlist.simplegraph_df <- function(x, ...) {

  adjlist <- modifyList(
    structure(
      replicate(nrow(x$nodes), character()),
      names = x$nodes[[1]]
    ),
    tapply(x$edges[,2], x$edges[,1], c, simplify = FALSE)
  )

  graph(adjlist)
}
