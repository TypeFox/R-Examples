#' Print for object of class PresenceAbsence
#' @author Bruno Vilela
#' 
#' @description Print for objects of class PresenceAbsence.
#' 
#' @usage 
#' \method{print}{PresenceAbsence}(x, \dots)
#' 
#' @param x an object of class \code{\link{PresenceAbsence}}.
#' @param ... Other print parameters.
#' 
#' @S3method print PresenceAbsence


print.PresenceAbsence <- function(x, ...) {
  resolution <- res(x$Ric)
  cat("\nClass:", class(x),
      "\nNumber of species:", (ncol(x$Pre) - 2),
      "\nNumber of cells:", nrow(x$Pre))
  cat("\nResolution: ", resolution[1], ", ", resolution[2], " (x, y)\n", sep="")  
  
}

