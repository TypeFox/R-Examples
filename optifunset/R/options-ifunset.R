#' Set Options if Unset
#' 
#' This function will set an option if it isn't already present within the global options 
#' returned by the \code{options()} function
#' 
#' @param ... any options can be defined, using \code{name = value}, if they are not already defined. 
#' 
#' Options can also be passed by giving a single unnamed argument which is a named list.
#' @param force Force the Option to Be Set
#' @examples
#' options.ifunset(width=10)            #IGNORED, ALREADY EXISTS
#' options.ifunset(width=10,force=TRUE) #FORCED UPDATE TO OPTION
#' options.ifunset(myoption=TRUE)       #New Option is Created
options.ifunset <- function(...,force=FALSE){
  d <- data.frame(..., stringsAsFactors = FALSE)
  o <- options()
  for(n in names(d)){ if(force | is.null(o[[n]])){ o[[n]] = d[1,n] }}
  options(o)
}