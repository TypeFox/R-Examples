#' @title Print gravity model
#' @description summary method for class "gravity"
#' @param x    Object of class gravity
#' @param ...  Ignored
#' @method print gravity
#' @export
"print.gravity" <- function(x, ...) { 
  cat("Gravity model\n\n") 
  print(summary(x$gravity))
  }
