#' WaterOneFlowNamespace
#'
#' A helper function that finds out the WaterOneFlow namespace information
#' based on the version number 1.0 or 1.1.
#'
#' @param version The version of the WaterOneFlow XML namespace. Must be either "1.0" or "1.1"
#' @return A list with the namespaces and corresponding prefixes. This namespace
#' information is important for correct parsing of the WaterML XML document.
#' @keywords WaterML
#' @export
#' @examples
#' ns <- WaterOneFlowNamespace("1.0")
#' ns <- WaterOneFlowNamespace("1.1")

WaterOneFlowNamespace <- function(version) {

  if (version == "1.0") {
    ns <- c(soap="http://schemas.xmlsoap.org/soap/envelope/",
            xsd="http://www.w3.org/2001/XMLSchema",
            xsi="http://www.w3.org/2001/XMLSchema-instance",
            sr="http://www.cuahsi.org/waterML/1.0/",
            gsr="http://www.cuahsi.org/his/1.0/ws/")
  } else {
    ns <- c(soap="http://schemas.xmlsoap.org/soap/envelope/",
            xsd="http://www.w3.org/2001/XMLSchema",
            xsi="http://www.w3.org/2001/XMLSchema-instance",
            sr="http://www.cuahsi.org/waterML/1.1/",
            gsr="http://www.cuahsi.org/his/1.1/ws/")
  }
  return(ns)
}
