#' Add basic metadata
#' 
#' adds Dublin Core metadata elements to (top-level) nexml 

#' @param title A title for the dataset
#' @param description a description of the dataset
#' @param creator name of the data creator.  Can be a string or R person object
#' @param pubdate publication date.  Default is current date.  
#' @param rights the intellectual property rights associated with the data.  
#'   The default is Creative Commons Zero (CC0) public domain declaration, 
#'   compatiable with all other licenses and appropriate for deposition 
#'   into the Dryad or figshare repositories.  CC0 is also recommended by the Panton Principles.
#'   Alternatively, any other plain text string can be added and will be provided as the content
#'   attribute to the dc:rights property.  
#' @param publisher the publisher of the dataset. Usually where a user may go to find the canonical
#'   copy of the dataset: could be a repository, journal, or academic institution.  
#' @param citation a citation associated with the data.  Usually an acompanying academic journal 
#'    article that indicates how the data should be cited in an academic context.  Multiple citations
#'    can be included here.  
#'    citation can be a plain text object, but is preferably an R `citation` or `bibentry` object (which
#'    can include multiple citations.  See examples
#' @param nexml a nexml object to which metadata should be added.  A new 
#'   nexml object will be created if none exists. 
#' @return an updated nexml object
#' @details \code{add_basic_meta()} is just a wrapper for \code{\link{add_meta}} to make it easy to 
#'    provide generic metadata without explicitly providing the namespace.  For instance, 
#'    \code{add_basic_meta(title="My title", description="a description")} is identical to:
#'    \code{add_meta(list(meta("dc:title", "My title"), meta("dc:description", "a description")))}

#'    Most function arguments are mapped directly to the Dublin Core terms
#'    of the same name, with the exception of `rights`, which by default maps
#'    to the Creative Commons namespace when using CC0 license.
#'    
#' @seealso \code{\link{add_trees}} \code{\link{add_characters}} \code{\link{add_meta}}
#' @export 
#' @importFrom stats na.omit
#' @importFrom utils capture.output head object.size
#' @examples
#' nex <- add_basic_meta(title = "My test title",
#'              description = "A description of my test",
#'              creator = "Carl Boettiger <cboettig@@gmail.com>",
#'              publisher = "unpublished data",
#'              pubdate = "2012-04-01")
#' 
#'  ## Adding citation to an R package:
#'  nexml <- add_basic_meta(citation=citation("ape"))
#' \dontrun{
#'  ## Use knitcitations package to add a citation by DOI:
#'  library(knitcitations)
#'  nexml <- add_basic_meta(citation=cite("10.2307/2408428"))
#'  }
#' @include classes.R
add_basic_meta <- function(title = NULL, 
                           description = NULL,
                           creator = Sys.getenv("USER"),
                           pubdate = Sys.Date(),
                           rights = "CC0",
                           publisher = NULL,
                           citation = NULL,
                           nexml = new("nexml")
                           ){

  mymeta <- get_metadata(nexml)
  m <- mymeta$content
  names(m) <- mymeta$property
  
  
  if(!is.null(title)) 
    nexml <- add_meta(meta("dc:title", title), nexml)
  if(!is.null(creator) || creator == "") 
    nexml <- add_meta(meta("dc:creator", format(creator)), nexml)
  if(!is.null(pubdate)) 
    if(!is.null(m)) 
      if(is.null(m["dc:pubdate"]) | is.na(m["dc:pubdate"]))
        nexml <- add_meta(meta("dc:pubdate", format(pubdate)), nexml)
  if(!is.null(description)) 
    nexml <- add_meta(meta("dc:description", description), nexml)
  if(!is.null(rights)){
    if(rights == "CC0")
      if(is.null(get_license(nexml)))
      nexml <- add_meta(meta(rel="cc:license", 
                             href="http://creativecommons.org/publicdomain/zero/1.0/"), nexml)
    else 
      nexml <- add_meta(meta("dc:rights", rights), nexml)
  }
  if(!is.null(citation))
    if(is(citation, "BibEntry"))
      class(citation) = "bibentry"
    if(is(citation, "bibentry"))
      nexml <- add_meta(nexml_citation(citation), nexml)
    else
      nexml <- add_meta(meta("dcterms:bibliographicCitation", citation), nexml)
  nexml
}

