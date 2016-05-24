
#' Breadth-first search of a graph
#'
#' @param graph Input graph.
#' @param from Character vector, which vertices to start the search
#'   from. By default all vertices are attempted.
#' @return Character vector of the named of the visited vertices,
#'   in the order of their visit.
#'
#' @export
#' @examples
#' funcs <- graph(list(
#'   drop_internal = character(0),
#'   get_deps = c("get_description", "parse_deps",
#'     "%||%", "drop_internal"),
#'   get_description = "pkg_from_filename",
#'   parse_deps = "str_trim",
#'   cran_file = c("get_pkg_type", "r_minor_version", "cran_file"),
#'   download_urls = c("split_pkg_names_versions", "cran_file"),
#'   filename_from_url = character(0),
#'   get_pkg_type = character(0),
#'   pkg_download = c("dir_exists", "download_urls",
#'     "filename_from_url", "try_download"),
#'   r_minor_version = character(0),
#'   try_download = character(0),
#'   drop_missing_deps = character(0),
#'   install_order = character(0),
#'   restore = c("pkg_download", "drop_missing_deps",
#'     "install_order", "get_deps"),
#'   snap = character(0),
#'   `%||%` = character(0),
#'   data_frame = character(0),
#'   dir_exists = character(0),
#'   pkg_from_filename = character(0),
#'   split_pkg_names_versions = "data_frame",
#'   str_trim = character(0)
#' ))
#' bfs(funcs)

bfs <- function(graph, from = vertex_ids(graph)) {

  graph <- as_graph_adjlist(graph)

  V <- names(graph)
  N <- length(V)
  result <- character()
  q <- character()
  marks <- structure(rep(FALSE, N), names = V)

  while (length(from)) {

    s <- from[1]
    from <- from[-1]
    if (!marks[[s]]) {
      result <- c(result, s)
      marks[[s]] <- TRUE
      q <- c(q, s)
    }

    while (length(q)) {

      s2 <- q[1]
      q <- q[-1]

      for (n in graph[[s2]]) {
        if (!marks[[n]]) {
          result <- c(result, n)
          q <- c(q, n)
          marks[[n]] <- TRUE
        }
      }
    }
  }

  result
}
