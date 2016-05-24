#' WaterOneFlowVersion
#'
#' A helper function that finds out the WaterOneFlow service version from
#' the URL of the wsdl file. By default it checks for "cuahsi_1_0" or "cuahsi_1_1" in
#' the url. If that is not found, then the function checks the version inside the
#' WSDL file.
#'
#' @param WSDL The URL of the WSDL, for example http://icewater.usu.edu/MudLake/cuahsi_1_0.asmx?WSDL
#' @return A list with two items: Version (either 1.0 or 1.1), and Namespace
#' (either http://www.cuahsi.org/his/1.0/ws/ or http://www.cuahsi.org/his/1.1/ws/)
#' @keywords WaterML
#' @export
#' @examples
#' versionInfo <- WaterOneFlowVersion("http://icewater.usu.edu/MudLake/cuahsi_1_0.asmx?WSDL")
#' versionInfo$Version
#' versionInfo$Namespace

WaterOneFlowVersion <- function(WSDL) {

  #check CUAHSINamespace parameter
  validNamespace_1_0 <- "http://www.cuahsi.org/his/1.0/ws/"
  validNamespace_1_1 <- "http://www.cuahsi.org/his/1.1/ws/"

  if (grepl("cuahsi_1_1", WSDL)) {
    return(list(Version="1.1", Namespace=validNamespace_1_1))
  }
  if (grepl("cuahsi_1_0", WSDL)) {
    return(list(Version="1.0", Namespace=validNamespace_1_0))
  }

  #if not found then assume 1.0
  return(list(Version="1.0", Namespace=validNamespace_1_0))
}
