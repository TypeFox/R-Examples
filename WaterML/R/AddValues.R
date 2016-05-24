#' AddValues
#'
#' This function adds a table of data values to HydroServer Lite.
#' The input must be a data.frame with Time and DataValue fields
#' The Time field must be POSIXct format and DataValue must be numeric format
#' it is also required to enter a valid SiteID, VariableID, SourceID, MethodID and
#' QualityControlLevelID. New data values shall be inserted only if the SiteID, VariableID,
#' SourceID, MethodID and QualityControlLevelID entries already exist in the HydroServer.
#' NOTE: this only works with HydroServer Lite that implements the JSON API.
#' you must specify a valid server url, user name, and password for the HydroServer.
#' The examples here use the 'sandbox' HydroServer on http://worldwater.byu.edu/app/
#' with the username: admin and password: password.
#'
#' @import RJSONIO
#' @import httr
#' @param server The URL of the web service ending with /services or with ?wsdl,
#'  for example: http://worldwater.byu.edu/app/index.php/default/services/cuahsi_1_1.asmx?wsdl
#'  alternatively you can specify the JSON API url like:
#'  http://worldwater.byu.edu/app/index.php/default/services/api/
#' @param username The valid HydroServer Lite username, for example "admin"
#' @param password The valid HydroServer Lite password, for example "password"
#' @param values The valid table of data values. This table must have the following columns:
#' Time (POSIXct), DataValue (numeric).
#' @param site The valid SiteID
#' @param variable The valid VariableID
#' @param methodID The valid MethodID
#' @param sourceID The valid SourceID
#' @param qualityControl The valid QualityControlLevelID
#' @return Status (the status showing if the values were added: OK or Error). If the status is Error, then
#' the Error message with reason why the values could not be added is also shown.
#' @keywords WaterML
#' @export
#' @examples
#' user <- "admin"
#' pass <- "password"
#' server <- "http://worldwater.byu.edu/app/index.php/default/services/cuahsi_1_1.asmx"
#' sourceID = 15
#' qualityID = 1
#' variableID = 43
#' siteID = 170
#' methodID = 10
#' random_times <- sort(Sys.time() + runif(3, 0, 10)*60)
#' random_values <- runif(3, 0, 100)
#' my_values <- data.frame(Time=random_times, DataValue=random_values)
#'
#' status  <- AddValues(server, username=user, password=pass,
#'                      site=siteID, variable=variableID,
#'                      methodID=methodID, sourceID=sourceID,
#'                      qualityControl=qualityID, values=my_values)

AddValues <- function(server, username, password, site, variable, methodID, sourceID,
                      qualityControl, values) {

  #check if the server is a valid url
  cuahsi <- regexpr("/cuahsi", server)
  services_api <- regexpr("/services/api", server)
  url <- NULL
  if (cuahsi > 0) {
    baseurl <- substr(server, 1, cuahsi)
    url <- paste(baseurl, "api/values",sep="")
  } else if (services_api > 0) {
    baseurl <- substr(server, 1, services_api)
    url <- paste(baseurl, "/services/api/values")
  } else {
    stop("The server url must contain cuahsi_1_1.asmx or ?wsdl or /services/api ")
  }

  #check if table has all required columns
  cols <- names(values)
  cols.required <- c("Time",
                     "DataValue"
                    )
  cols.matched <- match(cols.required, cols)
  if (length(cols.required[is.na(cols.matched)]) > 0) {
    cols.missing <- cols.required[is.na(cols.matched)]
    msg <- paste("values table has missing columns:", cols.missing)
    stop(msg)
  }

  i <- 1
  N <- nrow(values)

  dates_formatted <- format(values$Time, "%Y-%m-%d %H:%M:%S")
  value_list <- as.matrix(cbind(dates_formatted, values$DataValue))
  colnames(value_list) <- NULL
  rownames(value_list) <- NULL

  x <- list(
    user = username,
    password = password,
    SiteID = site,
    VariableID = variable,
    MethodID = methodID,
    SourceID = sourceID,
    values = value_list
  )
  post.body <- RJSONIO::toJSON(x)
  print(post.body)
  #post data to server:
  response <- POST(url,
                   body = post.body,
                   add_headers("Content-Type" = "application/json")
  )
  status.code <- http_status(response)$category
  if (tolower(status.code) == "success") {
    status = content(response, type="application/json")
    print(status)
    return(paste(status$status, status$message))
  } else {
      print(response)
      status = print(response)
  }
  return (status)
}
