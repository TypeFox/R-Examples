
#' Topological sorting of a graph
#'
#' @param graph Input graph.
#' @return Character vector of vertex ids, in topological order.
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
#' topological_sort(remove_loops(funcs))

topological_sort <- function(graph) {

  graph <- as_graph_adjlist(graph)

  V <- names(graph)
  N <- length(V)

  ## some easy cases
  if (length(graph) <= 1 ||
      sum(sapply(graph, length)) == 0) return(V)

  marked <- 1L; temp_marked <- 2L; unmarked <- 3L
  marks <- structure(rep(unmarked, N), names = V)
  result <- character(N)
  result_ptr <- N

  visit <- function(n) {
    if (marks[n] == temp_marked) stop("Call graph not a DAG: ", n)
    if (marks[n] == unmarked) {
      marks[n] <<- temp_marked
      for (m in graph[[n]]) visit(m)
      marks[n] <<- marked
      result[result_ptr] <<- n
      result_ptr <<- result_ptr - 1
    }
  }

  while (any(marks == unmarked)) {
    visit(names(which(marks == unmarked))[1])
  }

  result
}
