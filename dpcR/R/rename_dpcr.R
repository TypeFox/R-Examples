#' Rename object
#' 
#' Renames objects of class \code{\linkS4class{adpcr}} or \code{\linkS4class{ddpcr}}.
#' 
#' @param x an \code{\linkS4class{adpcr}} or \code{\linkS4class{ddpcr}} object.
#' @param exper a \code{factor} of new experiments' names. If \code{NULL},
#' experiments' names are not changed.
#' @param replicate a \code{factor} of new replicates' ids. If \code{NULL},
#' replicates' names are not changed.
#' @keywords manip
#' @export


rename_dpcr <- function(x, exper = NULL, replicate = NULL) {
  #add check if numeric
  if (!is.null(exper)) {
    if(length(exper) == 1) {
      exper <- as.factor(rep(exper, ncol(slot(x, ".Data"))))
    }
    slot(x, "exper") <- exper
  }
    
  if (!is.null(replicate))
    slot(x, "replicate") <- replicate
  
  colnames(slot(x, ".Data")) <- paste0(slot(x, "exper"), ".", slot(x, "replicate"))
  x
}