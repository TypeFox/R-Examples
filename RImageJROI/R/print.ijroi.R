##' @title Print \code{ijroi} objects
##' @param all logical indicating whether to print all information
##' from ijroi object as opposed to a subset of relevant information.
##' Defaults to \code{FALSE}.
##' @param x \code{ijroi} object to be printed.
##' @param ... further arguments passed to \code{\link{print}}.
##' @method print ijroi
##' @export
##' @author Mikko Vihtakari, David Sterratt
##' @seealso \code{\link{read.ijroi}}
print.ijroi <- function(x, all=FALSE, ...) {
  if (!all) {
    ## Exclude fields irrelevant to the type
    if (x$type %in% x$types[c("rect", "oval", "point")]) {
      exclude<- c("x1", "y1", "x2", "y2")
    } else {
      exclude <- c("bottom", "left", "top", "right", "width", "height")
    }
    x <- x[!names(x) %in% exclude]

    ## Exclude irrespective of the type
    x <- x[!names(x) %in% c("version", "types")]

    ## Exclude elements that are equal to 0
    exclude.if.0 <- c("n", "strokeWidth", "shapeRoiSize", "strokeColor",
                      "fillColor", "style", "headSize", "arcSize", "position")
    exclude.these <- unlist(lapply(x[names(x) %in% exclude.if.0],
                                  function(k) { c(k == 0 | is.na(k)) }))
    exclude.these <- names(exclude.these[exclude.these == TRUE])
    x <- x[!names(x) %in% exclude.these]
  }
  print(x, ...)
}
