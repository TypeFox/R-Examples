#' Print
#'
#' Prints geozoo objects with tourr or prints them
#'
#' @export
#' @method print geozoo
#' @param x geozoo object
#' @param ... other arguements
#' @references \url{http://schloerke.github.io/geozoo/}
#' @author Barret Schloerke
#' @examples
#' \dontrun{
#' example(boy.surface)
#' example(cube.face)
#' }
#' @keywords dynamic
print.geozoo <- function(x, ...){
  colnames(x$points) <- paste("V", seq_len(ncol(x$points)))


  if ( (! requireNamespace("tourr")) || (! interactive()) ) {
    as.list(x)
  } else {
    if (nrow(x$points) > 1000) {
      message("Reduce points to improve performance\n")
    }

    if (!is.null(x$edges)) {
      message("Edges will not be printed using 'tourr'\n")
    }

    tourr::animate(
      x$points,
      tour_path = tourr::grand_tour(),
      tourr::display_xy(...),
      ...
    )
  }

}



#' Print Without Rescale
#'
#' Prints objects without rescaling them to 0, 1 in each dim
#'
#' @export
#' @method print geozooNoScale
#' @param x geozoo object
#' @param ... other arguements
#' @references \url{http://schloerke.github.io/geozoo/}
#' @author Barret Schloerke
#' @examples
#' \dontrun{
#' torus()
#' }
#' @keywords dynamic
print.geozooNoScale <- function(x, ...) {
  class(x) <- class(x)[-1]
  print(x, rescale = FALSE, ...)
}
