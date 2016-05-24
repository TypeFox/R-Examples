#' @export
#' 
SDcheck.default <- function(keeper, string, summarize = F){
  stop("object \"",substitute(keeper),"\" is not of the appropriate class")
}