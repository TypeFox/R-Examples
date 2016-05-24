#' Opens a web page associated with a DOI (digital object identifier).
#'
#' Uses the DOI name of a study reference to locate the e-journal website,
#'    or reference/citation website in Web of Science or Google Scholar.  Opens
#'    in default web-browser.    
#'
#' @param theDOI A string that identifies an electronic document on the web.
#' @param host A string that defines the domain link used to open the DOI.  The 
#'    default, \code{"DOI"}, will open to the web page associated with 
#'    the DOI (e.g., publisher website).  Other options include \code{"WOS"} that 
#'    will open the DOI in Web of Science, and \code{"GS"} which will open the browser
#'    to Google Scholar.
#'
#' @return NULL
#'
#' @examples \dontrun{
#'
#' browse_DOI("10.1086/603628")        
#'}
#'
#' @importFrom utils browseURL
#' @export browse_DOI

browse_DOI <- function(theDOI, 
                       host = "DOI") {
                      
  theDomain <- switch(host, 
                      DOI = "http://dx.doi.org/",
                      WOS = "http://ws.isiknowledge.com/cps/openurl/service?url_ver=Z39.88-2004&rft_id=info:doi/",
                      GS  = "http://scholar.google.com/scholar?q="
  )
  browseURL(paste0(theDomain, theDOI))
}