#' Summary for object of class PresenceAbsence
#' @author Bruno Vilela
#' 
#' @description Summary for objects of class PresenceAbsence.
#' 
#' @usage 
#' \method{summary}{PresenceAbsence}(object, \dots)
#' 
#' @param object an object of class \code{\link{PresenceAbsence}}.
#' @param ... Other summary parameters.
#' 
#' @S3method summary PresenceAbsence

summary.PresenceAbsence <- function(object, ...) {
  
  class <- class(object)
  Numberofspecies <- ncol(object$Pre) - 2
  Numberofcells <- nrow(object$Pre)
  x2 <- object$Pre[, -(1:2), drop = FALSE]
  
  if (is.vector(x2)) {
    nomes <- names(x2)
    x2 <- matrix(x2, ncol = length(x2))
    colnames(x2) <- nomes          
  }
  
  Cellswithpresence <- sum(rowSums(x2) > 0)
  Cellswithoutanypresence <- sum(rowSums(x2) == 0)
  Specieswithoutanypresence <- sum(colSums(x2) == 0)
  SpeciesLargestRange <- names(2 + which(colSums(x2) == max(colSums(x2))))
  
  resolution <- res(object$Ri)
  extention <- extent(object$Ri)
  coordRef <- projection(object$R)  
  result <- list(class                      = class,
                 Numberofspecies            = Numberofspecies,
                 Numberofcells              = Numberofcells, 
                 Cellswithpresence          = Cellswithpresence, 
                 Cellswithoutanypresence    = Cellswithoutanypresence,
                 Specieswithoutanypresence  = Specieswithoutanypresence,
                 SpeciesLargestRange        = SpeciesLargestRange,
                 resolution                 = resolution,
                 extention                  = extention,
                 coordRef                   = coordRef)
  class(result) <- "summary.PresenceAbsence" 
  return(result)
}



