
#' add namespaces
#' 
#' add namespaces, avoiding duplication if prefix is already defined
#' @param namespaces a named character vector of namespaces
#' @param nexml a nexml object. will create a new one if none is given.  
#' @return a nexml object with updated namespaces 
#' @examples
#' ## Create a new nexml object with a single metadata element: 
#' modified <- meta(property = "prism:modificationDate", content = "2013-10-04")
#' nex <- add_meta(modified) # Note: 'prism' is defined in nexml_namespaces by default.  
#' 
#' ## Write multiple metadata elements, including a new namespace:  
#' website <- meta(href = "http://carlboettiger.info", 
#'                 rel = "foaf:homepage")              # meta can be link-style metadata
#' nex <- add_meta(list(modified,  website), 
#'                 namespaces = c(foaf = "http://xmlns.com/foaf/0.1/"))
#'
#' ## Append more metadata, and specify a level: 
#' history <- meta(property = "skos:historyNote",
#'                  content = "Mapped from the bird.orders data in the ape package using RNeXML")
#' nex <- add_meta(history, 
#'                 nexml = nex,
#'                 level = "trees",
#'                 namespaces = c(skos = "http://www.w3.org/2004/02/skos/core#"))
#' 
#' @seealso \code{\link{meta}} \code{\link{add_meta}}
#' @export 
add_namespaces <- function(namespaces, nexml = new("nexml")){
  if(!is.null(namespaces)){
## check for duplicated abbreviation, not for duplicated URI. OKAY to have multiple abbrs for same URI...

## FIXME Make sure that cases where abbreviation match actually match the URI as well
    notdups <- match(names(namespaces), names(nexml@namespaces)) 
    notdups <- sapply(notdups, is.na)
    if(all(notdups)) # all are unique 
      nexml@namespaces <-  c(nexml@namespaces, namespaces)
    else {
      nexml@namespaces <-  c(nexml@namespaces, namespaces[notdups])
    } 
  }
  nexml
}
