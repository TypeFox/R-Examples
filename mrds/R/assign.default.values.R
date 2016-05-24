#' Assign default values to list elements that have not been already assigned
#'
#' Assigns default values for \code{argument} in list \code{x} from
#' \code{argument=value} pairs in \dots{} if \code{x$argument} doesn't already
#' exist
#'
#' @param x generic list
#' @param \dots unspecified list of argument=value pairs that are used to
#'   assign values
#' @return x - list with filled values
#' @author Jeff Laake
#' @keywords ~utility
assign.default.values <- function(x,...){
  args <- list(...)
  notset <- !names(args)%in%names(x)
  x <- c(x,args[notset])
  return(x)
}
