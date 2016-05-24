#' originr - Species Origin Data
#'
#' @importFrom stats na.omit setNames
#' @importFrom httr GET content stop_for_status warn_for_status
#' @importFrom jsonlite fromJSON
#' @importFrom taxize get_uid classification get_tsn itis_native
#' @importFrom xml2 read_xml xml_find_all xml_text
#' @importFrom data.table rbindlist setDF
#' @name originr-package
#' @aliases originr
#' @docType package
#' @keywords package
#' @author Scott Chamberlain \email{myrmecocystus@@gmail.com}
#' @author Ignasi Bartomeus \email{nacho.bartomeus@@gmail.com}
#'
#' @section Data sources in the package:
#' \itemize{
#'  \item Encyclopedia of Life (http://eol.org)
#'  \item Flora Europaea (http://rbg-web2.rbge.org.uk/FE/fe.html)
#'  \item Global Invasive Species Database (http://www.iucngisd.org/gisd)
#'  \item Native Species Resolver (http://bien.nceas.ucsb.edu/bien/tools/nsr/nsr-ws/)
#'  \item Integrated Taxonomic Information Service (http://www.itis.gov/)
#' }
NULL
