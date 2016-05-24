#' R interface to the Biodiversity Heritage Library API.
#'
#' You need an API key to use the Biodiversity Heritage Library API. Get your
#' BHL API key at \url{http://www.biodiversitylibrary.org/getapikey.aspx}.
#' Put your API key in your .Rprofile file using e.g.,
#' `options(BioHerLibKey = "YOURBHLAPIKEY")`, and the functions within this package
#' will be able to use your API key without you having to enter it every time
#' you run a search.
#'
#' @importFrom httr GET content stop_for_status
#' @importFrom jsonlite fromJSON
#' @importFrom plyr ldply rbind.fill
#' @importFrom XML xmlSize xpathSApply xpathApply xmlParse xmlValue
#' @importFrom stats setNames
#' @importFrom methods is
#' @importFrom utils head
#' @name rbhl-package
#' @aliases rbhl
#' @docType package
#' @title R interface to the Biodiversity Heritage Library API.
#' @author Scott Chamberlain \email{myrmecocystus@@gmail.com}
#' @keywords package
NULL

#' Data.frame of all the BHL API methods from the BHL website.
#' @name rbhlmethods
#' @docType data
#' @keywords datasets
NULL
