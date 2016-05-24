#' NOT IN 
#' 
#' Commonly created NOT-IN operator
#' 
#' @param x object on the lhs 
#' @param table object/list on the rhs
#' 
#' @name notin
#' @rdname not-in
#' @export 

`%!in%` <- function(x,table) !`%in%`(x,table)
