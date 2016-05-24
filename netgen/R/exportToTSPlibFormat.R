#' Exports a network to the TSPlib format.
#'
#' @note Currently we only support euclidean 2D instances. Furthermore note, that
#' if \code{use.extended.format} is \code{TRUE}, most alternative TSPlib parsers
#' will most probably not be able to parse the generated file.
#'
#' @param x [\code{Network}]\cr
#'   Network to export.
#' @param filename [\code{character(1)}]\cr
#'   File name.
#' @param name [\code{character(1)} | \code{NULL}]\cr
#'   Character string describing the instance. Used for the NAME field in the
#'   TSPlib file format. Otherwise, the name of the instance is used. If the
#'   latter is \code{NULL}, this parameter is mandatory.
#' @param comment [\code{character(1)} | \code{NULL}]\cr
#'   Optional string with additional information about the instance. Used for
#'   the COMMENT field. If not provided the comment field of the instance is
#'   used. If the latter is \code{NULL}, no comment at all is saved.
#' @param use.extended.format [\code{logical(1)}]\cr
#'   Use the \dQuote{extended tsplib format} with additional information like cluster
#'   membership and bounds? Default is \code{TRUE}.
#' @param full.matrix [\code{logical(1)}]\cr
#'   Make use of \dQuote{FULL\_MATRIX} \dQuote{EDGE\_WEIGHT\_FORMAT} instead of
#'   node coordinates?
#'   Default is \code{FALSE}.
#' @param digits [\code{integer(1)}]\cr
#'   Round coordinates to this number of digits. Default is 10.
#' @return Nothing
#' @export
exportToTSPlibFormat = function(x, filename,
  name = NULL, comment = NULL,
  use.extended.format = TRUE,
  full.matrix = FALSE,
  digits = 10L) {
  assertFlag(full.matrix)
  if (is.null(name) && is.null(x$name)) {
    stopf("Please provide a name for the instance via the 'name' parameter.")
  }
  if (hasDepots(x)) {
    stopf("Currently only instances without depots can be exported to the tsplib format.")
  }
  name = BBmisc::coalesce(name, x$name)
  comment = BBmisc::coalesce(comment, x$comment)
  coordinates = x$coordinates
  n = nrow(coordinates)
  n.cluster = getNumberOfClusters(x)
  out = paste0("NAME : ", name, "\n")
  if (!is.null(comment)) {
    for (com in comment) {
      out = paste0(out, "COMMENT : ", com, "\n")
    }
  }
  out = paste0(out, "TYPE : TSP\n")
  out = paste0(out, "DIMENSION : ", n, "\n")
  if (!full.matrix) {
    out = paste0(out, "EDGE_WEIGHT_TYPE : ", if (is.null(x$edge.weight.type)) "EUC_2D" else x$edge.weight.type, "\n")
    if (use.extended.format) {
      out = paste0(out, "LOWER : ", x$lower, "\n")
      out = paste0(out, "UPPER : ", x$upper, "\n")
    }
    out = paste0(out, "NODE_COORD_SECTION\n")
    #FIXME: this works only for the 2d case
    for (i in seq(n)) {
      out = paste0(out, i, " ",
        round(coordinates[i, 1], digits = digits), " ",
        round(coordinates[i, 2], digits = digits),
        if (i < n || n.cluster > 1L) "\n" else "")
    }
    if (use.extended.format & n.cluster > 1L) {
      out = paste0(out, "CLUSTER_MEMBERSHIP_SECTION\n")
      membership = x$membership
      for (i in seq(n)) {
        out = paste0(out, membership[i], if (i < n) "\n" else "")
      }
    }
  } else {
    out = paste0(out, "EDGE_WEIGHT_TYPE : EXPLICIT\n")
    out = paste0(out, "EDGE_WEIGHT_FORMAT : FULL_MATRIX\n")
    out = paste0(out, "EDGE_WEIGHT_SECTION\n")
    for (i in seq(n)) {
      out = paste(out, collapse(x$distance.matrix[i, ], sep = " "), if (i < n) "\n" else "")
    }
  }
  out = paste(out, "\nEOF", sep = "")
  write(x = out, file = filename)
}
