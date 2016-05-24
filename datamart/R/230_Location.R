#' S4 base class to represent output Locations
#'
#' The \code{Location} class is an abstract class that 
#' represents a place where to put something.
#'
#' @name Location-class
#' @rdname Location-class
#' @exportClass Location
setClass(Class="Location", representation=representation(), contains="Xdata")

#' @rdname queries-methods
#' @name queries
#' @export
#' @docType methods
#' @aliases queries queries,Location-method
setMethod(
  f="queries",
  signature=c("Location"),
  definition=function(self) c(rownames(meta(self)), callNextMethod())
)

#' S4 base class to represent output in Memory
#'
#' The \code{MemoryLocation} class represents the place in the memory of
#' the current R process.
#' 
#' @examples
#' getSlots("MemoryLocation")
#'
#' @name MemoryLocation-class
#' @rdname MemoryLocation-class
#' @exportClass MemoryLocation
setClass(Class="MemoryLocation", representation=representation(), contains="Location")

