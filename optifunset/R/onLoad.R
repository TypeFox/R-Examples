#' OptIfUnset Package
#' 
#' @description Set an option only if is not currently set
#' @details
#' This package contains a single function which only updates the global options if the option is presently unset
#' @examples
#' options.ifunset(width=100)  #NO CHANGE, ALREADY EXISTS
#' options.ifunset(myuniqueoption=TRUE) #NEW Option Created
#' @name optifunset
NULL