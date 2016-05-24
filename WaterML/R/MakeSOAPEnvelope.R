#' MakeSOAPEnvelope
#'
#' A helper function that makes a SOAP envelope to send to the
#' CUAHSI WaterOneFlow SOAP web service. It is internally used by the GetSites,
#' GetSiteInfo, GetVariables and GetValues functions.
#'
#' @import httr
#' @param CUAHSINamespace The SOAP namespace. This must be either "http://www.cuahsi.org/his/1.0/ws"
#' for WaterML 1.0, or "http://www.cuahsi.org/his/1.1/ws" for WaterML 1.1
#' @param MethodName The name of the WaterOneFlow web service method. It can be one of the following
#' values: "GetSites", "GetSitesObject", "GetSitesByBoxObject", "GetSiteInfoObject",
#' "GetVariablesObject", "GetValuesObject"
#' @param parameters An optional vector of named parameters for the web method. For GetSites,
#' GetSitesObject and GetVariables no parameters are required. For GetSiteInfoObject you need the
#' "site" parameter. For GetValuesObject you need the "location", "variable", "startDate" and "endDate"
#' parameters.
#' @return A <soap:Envelope> text in XML format. This text is send in a HTTP POST body to the
#' SOAP service. Two headers must be sent in the request: Content-Type="text/XML" and
#' SOAPAction=paste(CUAHSINamespace, MethodName). For example if MethodName is GetSites and
#' the WaterML version is 1.1, then SOAPAction="http://www.cuahsi.org/his/1.1/ws/GetSites".
#' @keywords WaterML
#' @export
#' @examples
#' library(httr)
#' library(XML)
#' myEnvelope <- MakeSOAPEnvelope("http://www.cuahsi.org/his/1.1/ws/", "GetSitesObject")
#' SOAPAction <- "http://www.cuahsi.org/his/1.1/ws/GetSitesObject"
#' url <- "http://hydrodata.info/chmi-d/cuahsi_1_1.asmx"
#' response <- POST(url, body = myEnvelope,
#'                  add_headers("Content-Type" = "text/xml", "SOAPAction" = SOAPAction),
#'                  verbose())
#' status.code <- http_status(response)$category
#' WaterML <- xmlParse(response)
#' WaterML

MakeSOAPEnvelope <- function(CUAHSINamespace, MethodName, parameters=NULL) {

  #check CUAHSINamespace parameter
  validNamespace_1_0 <- "http://www.cuahsi.org/his/1.0/ws/"
  validNamespace_1_1 <- "http://www.cuahsi.org/his/1.1/ws/"
  if (CUAHSINamespace != validNamespace_1_0 & CUAHSINamespace != validNamespace_1_1) {
    stop(paste("The CUAHSINamespace parameter must be either",
               validNamespace_1_0, "or", validNamespace_1_1))
  }
  #check MethodName
  validMethodNames <- c("GetSites", "GetSitesObject", "GetSitesByBoxObject",
                        "GetSiteInfoObject", "GetVariableInfoObject", "GetValuesObject")
  if (length(validMethodNames[validMethodNames == MethodName]) == 0) {
    message <- paste(validMethodNames, collapse=", ")
    stop(paste("The MethodName must be one of the following:", message))
  }

  soapAction <- paste(CUAHSINamespace, MethodName, sep="")
  #make the XML for parameters
  if (MethodName == "GetSitesObject" | MethodName == "GetSites") {
    XMLParameters <- c('<site></site>')
  } else {
    XMLParameters <- paste('<',names(parameters), '>', parameters, '</', names(parameters),'>\n',sep="")
  }
  envelope <- paste(
    '<?xml version="1.0" encoding="utf-8"?>\n',
    '<soap:Envelope xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">\n',
    '<soap:Body>\n',
    '<', MethodName, ' xmlns="', CUAHSINamespace, '">\n',
    paste(XMLParameters, sep="", collapse=""),
    '<authToken></authToken>\n',
    '</', MethodName, '>\n',
    '</soap:Body>\n',
    '</soap:Envelope>',
    sep="")

  return(envelope)
}
