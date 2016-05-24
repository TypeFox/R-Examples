#' oai
#'
#' @section Introduction:
#' oai is an R client to work with OAI-PMH (Open Archives Initiative Protocol
#' for Metadata Harvesting) services, a protocol developed by the
#' [Open Archives Initiative](https://en.wikipedia.org/wiki/Open_Archives_Initiative).
#' OAI-PMH uses XML data format transported over HTTP.
#'
#' @section OAI-PMH Info:
#' See the OAI-PMH V2 specification at
#' \url{http://www.openarchives.org/OAI/openarchivesprotocol.html}
#'
#' @section oai package details:
#' oai is built on \code{xml2} and `httr`. In addition, we give back data.frame's
#' whenever possible to make data comprehension, manipulation, and visualization
#' easier. We also have functions to fetch a large directory of OAI-PMH services -
#' it isn't exhaustive, but does contain a lot.
#'
#' @section Paging:
#' Instead of paging with e.g., \code{page} and \code{per_page} parameters, OAI-PMH uses
#' (optionally) \code{resumptionTokens}, with an optional expiration date. These tokens
#' can be used to continue on to the next chunk of data, if the first request did not
#' get to the end. Often, OAI-PMH services limit each request to 50 records, but this
#' may vary by provider, I don't know for sure. The API of this package is such that
#' we \code{while} loop for you internally until we get all records. We may in the future
#' expose e.g., a \code{limit} parameter so you can say how many records you want, but we
#' haven't done this yet.
#'
#' @name oai-package
#' @aliases oai
#' @importFrom methods is
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom httr GET content stop_for_status
#' @importFrom xml2 read_xml xml_children xml_text as_list xml_attrs xml_name xml_attr
#' @importFrom plyr rbind.fill
#' @importFrom stringr str_extract
#' @docType package
#' @title OAI-PMH Client
#' @author Scott Chamberlain \email{myrmecocystus@@gmail.com}
#' @keywords package
NULL

#' Metadata providers data.frame.
#'
#' @name providers
#' @docType data
#' @keywords datasets
#' @return A data.frame of three columns:
#' \itemize{
#'  \item repo_name - Name of the OAI repository
#'  \item base_url - Base URL of the OAI repository
#'  \item oai_identifier - OAI identifier for the OAI repository
#' }
NULL
