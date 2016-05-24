
#' Simple Graph Data Types and Basic Algorithms
#'
#' Simple classic graph algorithms for simple graph classes.
#' Graphs may possess vertex and edge attributes. 'simplegraph' has
#' no dependencies and it is writting entirely in R, so it is easy to
#' install.
#'
#' @docType package
#' @name simplegraph
NULL

#' Create a graph
#'
#' Graphs can be specified as adjacency lists or (two) data frames.
#'
#' If the first argument is a data frame, then it is interpreted as
#' vertex data, and a second data frame must be supplied as edge data.
#' The first column of the vertex data must contain (character) vertex
#' ids. The first two columns of the edge data frame must contain the
#' directed edges of the graph, in the order of tail and head, as
#' characters referring to the nodes ids. Other columns are kept as
#' metadata.
#'
#' If the first argument is not a data frame, but a list, then it is
#' interpreted as an adjacency list. It must be named, and the names
#' will be used as vertex ids. Each list element must be a character
#' vector containing the successors of each vertex.
#'
#' @param x A data frame, or a named list of character vectors. See
#'   details below.
#' @param ... Additional arguments, see details below.
#' @return A graph object.
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
#' funcs
#'
#' vertices <- data.frame(
#'   stringsAsFactors = FALSE,
#'   name = c("Tom Hanks", "Cate Blanchett", "Matt Damon", "Kate Winslet",
#'     "Saving Private Ryan", "Contagion", "The Talented Mr. Ripley"),
#'   what = c("actor", "actor", "actor", "actor", "movie", "movie", "movie"),
#'   born = c("1956-07-09", "1966-05-26", "1970-10-08", "1975-10-05",
#'     NA, NA, NA),
#'   gender = c("M", "F", "M", "F", NA, NA, NA),
#'   year = c(NA, NA, NA, NA, 1998, 2011, 1999)
#' )
#'
#' edges <- data.frame(
#'   stringsAsFactors = FALSE,
#'   actor = c("Tom Hanks", "Cate Blanchett", "Matt Damon", "Matt Damon",
#'     "Kate Winslet"),
#'   movie = c("Saving Private Ryan", "The Talented Mr. Ripley",
#'     "Saving Private Ryan", "The Talented Mr. Ripley", "Contagion")
#' )
#' actors <- graph(vertices, edges)
#' actors

graph <- function(x, ...)
  UseMethod("graph")

#' Check the validity of a graph data structure
#'
#' This is mainly for internal checks, but occasionally it
#' might also useful externally.
#'
#' @param x Graph.
#' @param ... Extra arguments are curently ignored.
#'
#' @export
#' @examples
#' G <- graph(list(A = c("B", "C"), B = "C", C = "A"))
#' sanitize(G)
#'
#' G <- c(G, list("this is not good" = c(1, 2, 3)))
#' try(sanitize(G))

sanitize <- function(x, ...)
  UseMethod("sanitize")
