#' Basic Mediawik API with one task to retrieve descriptions
#'
#' This function queries Wikipedia in the language of your choice. You
#' can select the language per two character country code Wikipedia uses in
#' they sub domains.
#'
#' @param term Is the term or string you are searching for on Wikipedia.
#' @param language Is the language of the Wikipedia (default: en). You
#'  	  can use two character country codes that Wikipedia uses in their
#' 	  sub domains.
#' @return The function returns a vector with description strings.
#' @examples \dontrun{
#' 		descriptions = bef.mediawiki.api.define(term = "Tree")
#' 		descriptions = bef.mediawiki.api.define(term = "Baum", language = "de")
#'	     }
#' @import RCurl
#' @import XML
#' @export bef.mediawiki.api.define

bef.mediawiki.api.define <- function(term, language="en") {
      wiki_api_url= paste0("http://", language, ".wikipedia.org/w/api.php")
	search_return <- getForm(
	  wiki_api_url,
	  action  = "opensearch",
	  search  = term,
	  format  = "xml",
	  .opts   = ""
	)
	document = xmlTreeParse(search_return, useInternalNodes=TRUE)
	document_root = xmlRoot(document)
	nodeset = getNodeSet(document_root, "//xmlns:Description", "xmlns")
	descriptions = xmlSApply(nodeset, xmlValue)
	cleaned_descriptions = sapply(descriptions, function(x) clean_html_string(x), USE.NAMES=F)
	return(cleaned_descriptions)
}
