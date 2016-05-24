#' WaterMLVersion
#'
#' A helper function that finds out the WaterML version from
#' the WaterML document. By default it checks for "http://www.opengis.net/waterml/2.0"
#' Otherwise it tries to detect "http://www.cuahsi.org/waterML/1.1/" (for WaterML 1.1) or
#' "http://www.cuahsi.org/WaterML/1.0/" (for WaterML 1.0)
#'
#' @param doc The XML document object
#' @return A character with the WaterML version: either 1.0, 1.1, or 2.0
#' @keywords WaterML
#' @export
#' @examples
#' library(httr)
#' library(XML)
#' url <- "http://www.waterml2.org/KiWIS-WML2-Example.wml"
#' response <- GET(url)
#' doc <- xmlParse(response)
#' version <- WaterMLVersion(doc)

WaterMLVersion <- function(doc) {

  #check namespaces
  ns_list <- xmlNamespaceDefinitions(doc, simplify=TRUE)
  namespaces <- as.character(unlist(ns_list))

  wml_2_0_namespace <- "http://www.opengis.net/waterml/2.0"
  wml_1_1_namespace <- "http://www.cuahsi.org/waterML/1.1/"
  wml_1_0_namespace <- "http://www.cuahsi.org/waterML/1.0/"

  if (wml_2_0_namespace %in% namespaces) {
    return ("2.0")
  }

  if (wml_1_1_namespace %in% namespaces) {
    return ("1.1")
  }

  if (wml_1_0_namespace %in% namespaces) {
    return ("1.0")
  }

  #if not found assume 1.1
  return ("1.1")
}
