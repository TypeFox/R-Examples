#' cut edges
#'
#' @param x sg graph object
#' @param data point pattern used for computing g
#' @param R cutting length
#' @param ... ignored
#'
#' Removes edges with length > R.
#'
#'
#' @export

cut.sg <- function(x, data, R, ...) {
  if(missing(R)) stop("Cutting length R must be given.")
  if(!is(x, "sg")) stop("g is not sg-object.")

  en <- cut_c(x$edges, sg_parse_coordinates(data), R)
  x$edges <- en
  x$note <- c(x$note, paste0("cut with R=", R))
  x
}
