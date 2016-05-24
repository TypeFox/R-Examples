#' Print without attributes
#' 
#' @param x The object to be printed out
#' @param list_attributes prints a list of attribute names if TRUE (default=TRUE)
#' @details To suppress the printout of attributes like "deriv". 
#' @export
print0 <- function(x, list_attributes=TRUE ) {
  attributes_all <- names(attributes(x))
  attributes_rm <- attributes_all[!(attributes_all %in% c("dim","names","dimnames","row.names","col.names"))]
  attributes(x)[attributes_rm] <- NULL
  print.default(x)
  if(list_attributes)
    cat("Attributes:",attributes_all)
}

