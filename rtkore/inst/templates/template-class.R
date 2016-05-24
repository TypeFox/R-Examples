#' The [\code{\linkS4class{EXAMPLE}}] class
#'
#' This class contains an example. This line goes into the description
#'
#' This line and the next ones go into the details.
#' This line thus appears in the details as well.
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{slot1}:}{Matrix of class \code{"numeric"}, containing data from slot1}
#'    \item{\code{slot2}:}{Object of class \code{"character"}, containing data that needs to go in slot2.}
#'  }
#'@examples
#' getSlots("EXAMPLE")
#' new("EXAMPLE")
#' 
#' @note You can still add notes
#' 
#' @author Serge Iovleff
#' 
#' @name EXAMPLE 
#' @rdname EXAMPLE
#' @aliases EXAMPLE-class
#' @exportClass EXAMPLE
#'
setClass("EXAMPLE",
    representation(slot1="numeric",
        slot2="character",
        myslot3="data.frame"),
    contains = "character"
)
