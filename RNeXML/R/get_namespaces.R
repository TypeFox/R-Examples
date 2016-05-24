
#' get namespaces
#' 
#' get namespaces
#' @param nexml a nexml object
#' @return a named character vector providing the URLs defining each
#' of the namespaces used in the nexml file.  Names correspond to 
#' the prefix abbreviations of the namespaces. 
#' @examples
#' comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
#' nex <- nexml_read(comp_analysis)
#' get_namespaces(nex)
#' @export 
get_namespaces <- function(nexml){
  nexml@namespaces
}
