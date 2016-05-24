#' @export
print.diffQ2object <- function(x, ...) {
  converted.x <- x
  class(converted.x) <- "diffQobject"
  print(converted.x)
  cat("Calculated 'left' Tm:", x[["xTm1.2.D2"]][1], "\n")
  cat("Calculated 'left' signal height:", x[["yTm1.2.D2"]][1], "\n")
  cat("Calculated 'right' Tm:", x[["xTm1.2.D2"]][2], "\n")
  cat("Calculated 'right' signal height:", x[["yTm1.2.D2"]][2], "\n")
}


