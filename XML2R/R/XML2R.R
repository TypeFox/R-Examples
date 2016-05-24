#' Parse XML files into (a list of) matrices or data frame(s)
#' 
#' This function is an experimental wrapper around \link{XML2Obs}. One should only use this function over \link{XML2Obs} if 
#' keys already exist in the XML data and ancestory doesn't need to be altered.
#' 
#' @param urls character vector or list of urls that point to an XML file (or anything readable by \link{xmlParse}).
#' @param xpath XML XPath expression that is passed to \link{getNodeSet}. If missing, the entire root and all descendents are captured and returned (ie, tables = "/"). 
#' @param df logical. Should matrices be coerced into data frames?
#' @return Returns list with one element for each relevant XML node. Each element contains a matrix by default.
#' @seealso \link{urlsToDocs}, \link{docsToNodes}, \link{nodesToList}, \link{listsToObs}
#' @export
#' @examples
#' \dontrun{
#' urls2 <- c("http://gd2.mlb.com/components/game/mlb/year_2013/mobile/346180.xml",
#'            "http://gd2.mlb.com/components/game/mlb/year_2013/mobile/346188.xml")
#' dat3 <- XML2R(urls2)
#' 
#' cens <- "http://www.census.gov/developers/data/sf1.xml"
#' census <- XML2R(cens)
#' }
#' 

XML2R <- function(urls, xpath, df=FALSE) {
  obs <- XML2Obs(urls, xpath)
  #add an option to try to automatically generate keys?
  collapse_obs(obs)
}

#' Parse XML files into a list of "observations"
#' 
#' This function takes a collection of urls that point to XML files and coerces the relevant info into a list of observations.
#' An "observation" is defined as a matrix with one row. An observation can also be thought of as a single instance of 
#' XML attributes (and value) for a particular level in the XML hierarchy. The names of the list reflect the XML node 
#' ancestory for which each observation was extracted from.
#' 
#' It's worth noting that a "url_key" column is appended to each observation to help track the origin of each observation.
#' The value of the "url_key" column is not the actual file name, but a simplified identifier to avoid unnecessarily repeating 
#' long file names for each observation. For this reason, an addition element (named "url_map") is added to the list of observations
#' in case the actual file named want to be used.
#' 
#' @param urls character vector or list of urls that point to an XML file (or anything readable by \link{xmlParse}).
#' @param xpath XML XPath expression that is passed to \link{getNodeSet}. If missing, the entire root and all descendents are captured and returned (ie, tables = "/"). 
#' @param append.value logical. Should the XML value be appended for relevant observations?
#' @param as.equiv logical. Should observations from two different files (but the same ancestory) have the same name returned?
#' @param url.map logical. If TRUE, the 'url_key' column will contain a condensed url identifier (for each observation)
#' and full urls will be stored in the "url_map" element. If FALSE, the full urls are included (for each observation) 
#' as a 'url' column and no "url_map" is included.
#' @param async logical. Allows for asynchronous download requests. This option is passed to the \code{async} option in the \code{RCurl::getURL} function.
#' @param quiet logical. Print file name currently being parsed?
#' @seealso \link{urlsToDocs}, \link{docsToNodes}, \link{nodesToList}, \link{listsToObs}
#' @return A list of "observations" and (possibly) the "url_map" element. 
#' @export
#' @examples
#' 
#' \dontrun{
#' urls <- c("http://gd2.mlb.com/components/game/mlb/year_2013/mobile/346180.xml",
#'            "http://gd2.mlb.com/components/game/mlb/year_2013/mobile/346188.xml")
#' obs <- XML2Obs(urls)
#' table(names(obs))
#' }

XML2Obs <- function(urls, xpath, append.value=TRUE, as.equiv=TRUE, url.map=FALSE, async=FALSE, quiet=FALSE) {
  if (missing(xpath)) xpath <- "/"  #select the root
  docs <- urlsToDocs(urls, async=async, quiet=quiet)
  valid.urls <- sapply(docs, function(x) attr(x, "XMLsource"))
  nodes <- docsToNodes(docs, xpath) 
  rm(docs)
  gc()
  l <- nodesToList(nodes)
  rm(nodes)
  gc()
  obs <- listsToObs(l, urls=valid.urls, append.value=append.value, as.equiv=as.equiv, url.map=url.map)
  return(obs)
}
