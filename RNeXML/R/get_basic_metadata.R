
## Goodness, but XPATH is so much more expressive for this purpose...
## get all top-level metadata  More extensible than hardwired functions
## The following methods are somewhat too rigid.  Might make more sense to do get_metadata(nexml, "nexml")["dc:creator"], etc.  
## Note that we define our namespace prefixes explicitly, so that should the NeXML use a different abberivation, this should still work.  

#' get_citation
#' 
#' get_citation
#' @param nexml a nexml object
#' @return the list of taxa
#' @export
get_citation <- function(nexml){
  b <- setxpath(as(nexml, "XMLInternalElementNode"))
## FIXME should return a citation class nexml! 
  cat(unname(xpathSApply(b, "/nexml/meta[@property='dcterms:bibliographicCitation']/@content", namespaces = nexml_namespaces)))
}

#' get_license
#'
#' get_license
#' @param nexml a nexml object
#' @return the list of taxa
#' @export 
get_license <- 
  function(nexml){
    b <- setxpath(as(nexml, "XMLInternalElementNode"))
    dc_rights <- unname(xpathSApply(b, "/nexml/meta[@property='dc:rights']/@content", namespaces = nexml_namespaces))
    cc_license <- unname(xpathSApply(b, "/nexml/meta[@rel='cc:license']/@href", namespaces = nexml_namespaces))
  if(length(dc_rights) > 0)
    dc_rights
  else
    cc_license
}






