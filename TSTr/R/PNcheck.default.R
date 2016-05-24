#' @export
#' 
PNcheck.default <- function(tree, string, useUpper = F){
  stop("object \"",substitute(tree),"\" is not of the appropriate class \"tstTree\"")
}