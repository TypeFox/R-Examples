#' @title API Authentication
#' 
#'  @description Initiate the Criteo API authentication process. The function returns the authentication token.
#'  
#'   @param user Username
#'   @param password Password
#'   @param company Company Name
#'   @param app App Name of API project
#'   @param version API version
#'   
#'   @export
#'   
#'   @return AuthToken
doCriteoAuth <- function(user, password, company, app, version){
  headerFields =
    c(Accept = "text/xml",
      Accept = "multipart/*",
      'Content-Type' = "text/xml; charset=utf-8",
      SOAPAction = "https://advertising.criteo.com/API/v201010/clientLogin")
  
  body = paste('<?xml version="1.0" encoding="utf-8"?>
               <soap:Envelope xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
               <soap:Body>
               <clientLogin xmlns="https://advertising.criteo.com/API/v201010">
               <username>', user ,'</username>
               <password>', password,'</password>
               <source>',company,'-', app ,'-', version,'</source>
               </clientLogin>
               </soap:Body>
               </soap:Envelope>',sep="")
  
  h = RCurl::basicTextGatherer()
  RCurl::curlPerform(url = "https://advertising.criteo.com/API/v201010/AdvertiserService.asmx",
              httpheader = headerFields,
              postfields = body,
              writefunction = h$update
  )
  
  authToken <- XML::xmlParse(h$value())
  authToken <- XML::xmlToList(authToken)
  authToken <- authToken$Body$clientLoginResponse$clientLoginResult
  
  authToken
}