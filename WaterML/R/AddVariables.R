#' AddVariables
#'
#' This function adds a table of variables to HydroServer Lite.
#' The input must be a data.frame with all required ODM 'variable' fields
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
#' @param variables The valid table of variables. This table must have the following columns:
#' VariableCode, VariableName, Speciation, VariableUnitsID, SampleMedium, ValueType,
#' IsRegular, TimeSupport, TimeUnitsID, DataType, GeneralCategory, NoDataValue.
#' NOTE that the values of these fields must be in the CUAHSI controlled vocabulary.
#' @return A table of the added variables, with two extra columns:
#' VariableID (the ID assigned by the server),
#' Status (the status showing if the variable was added: OK or Error). If the status is Error, then
#' the Error message with reason why the variable could not be added is also shown.
#' @keywords waterml
#' @export
#' @examples
#' user <- "admin"
#' pass <- "password"
#' server <- "http://worldwater.byu.edu/app/index.php/default/services/cuahsi_1_1.asmx"
#' #make random site codes
#' random_code = sprintf("R-%04d",sample(1:10000, 1))
#' my_variables <- data.frame(
#'   VariableCode = random_code,
#'   VariableName = "Color",
#'   Speciation = "Not Applicable",
#'   VariableUnitsID = 189,
#'   SampleMedium = "Groundwater",
#'   ValueType = "Sample",
#'   IsRegular = 1,
#'   TimeSupport = 0,
#'   TimeUnitsID = 100,
#'   DataType = "Average",
#'   GeneralCategory = "Hydrology",
#'   NoDataValue = -9999
#' )
#'
#' added_variables <- AddVariables(server, username=user, password=pass,
#'                                 variables=my_variables)

AddVariables <- function(server, username, password, variables) {

  #check if the server is a valid url
  cuahsi <- regexpr("/cuahsi", server)
  services_api <- regexpr("/services/api", server)
  url <- NULL
  if (cuahsi > 0) {
    baseurl <- substr(server, 1, cuahsi)
    url <- paste(baseurl, "api/variables",sep="")
  } else if (services_api > 0) {
    baseurl <- substr(server, 1, services_api)
    url <- paste(baseurl, "/services/api/variables")
  } else {
    stop("The server url must contain cuahsi_1_1.asmx or ?wsdl or /services/api ")
  }

  #check if variables table has all required columns
  cols <- names(variables)
  cols.required <- c("VariableCode",
                     "VariableName",
                     "Speciation",
                     "VariableUnitsID",
                     "SampleMedium",
                     "ValueType",
                     "IsRegular",
                     "TimeSupport",
                     "TimeUnitsID",
                     "DataType",
                     "GeneralCategory",
                     "NoDataValue")
  cols.matched <- match(cols.required, cols)
  if (length(cols.required[is.na(cols.matched)]) > 0) {
    cols.missing <- cols.required[is.na(cols.matched)]
    msg <- paste("variables table has missing columns:", cols.missing)
    stop(msg)
  }

  i <- 1
  N <- nrow(variables)
  added.variables <- data.frame(VariableID=character(),
                            VariableName=character(),
                            VariableCode=character(),
                            Speciation=character(),
                            VariableUnitsID=numeric(),
                            SampleMedium=character(),
                            ValueType=character(),
                            IsRegular=character(),
                            TimeSupport=numeric(),
                            TimeUnitsID=numeric(),
                            DataType=character(),
                            GeneralCategory=character(),
                            NoDataValue=numeric(),
                            Status=character(),
                            stringsAsFactors=FALSE)

  for(i in 1:N) {
    x <- list(
      user = username,
      password = password,
      VariableCode = as.character(variables$VariableCode[i]),
      VariableName = as.character(variables$VariableName[i]),
      Speciation=as.character(variables$Speciation[i]),
      VariableUnitsID=as.numeric(variables$VariableUnitsID[i]),
      SampleMedium=as.character(variables$SampleMedium[i]),
      ValueType=as.character(variables$ValueType[i]),
      IsRegular=as.character(variables$IsRegular[i]),
      TimeSupport=as.numeric(variables$TimeSupport[i]),
      TimeUnitsID=as.numeric(variables$TimeUnitsID[i]),
      DataType=as.character(variables$DataType[i]),
      GeneralCategory=as.character(variables$GeneralCategory[i]),
      NoDataValue=as.numeric(variables$NoDataValue[i]),
      Status=character()
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
      response_status = content(response, type="application/json")
      print(response_status)
      if (response_status$status == "200 OK") {
        new.id.start <- regexpr("ID=", response_status$message)
        new.id <- substr(response_status$message, new.id.start + 3, nchar(response_status$message))
        new.variable <-
          c(new.id,
            as.character(x$VariableName[i]),
            as.character(x$VariableCode[i]),
            as.character(x$Speciation[i]),
            as.numeric(x$VariableUnitsID[i]),
            as.character(x$SampleMedium[i]),
            as.character(x$ValueType[i]),
            as.character(x$IsRegular[i]),
            as.numeric(x$TimeSupport[i]),
            as.numeric(x$TimeUnitsID[i]),
            as.character(x$DataType[i]),
            as.character(x$GeneralCategory[i]),
            as.numeric(x$NoDataValue[i]),
            as.character(status.code)
        )
        added.variables[nrow(added.variables)+1,] <- new.variable
      } else {
        print(paste("ERROR!", status.code))
      }
      #http status code is other than success...
    } else {
      if (status.code == "server error") {
        print(response)
      }
      new.variable <-
        c(NA,
          as.character(x$VariableName[i]),
          as.character(x$VariableCode[i]),
          as.character(x$Speciation[i]),
          as.numeric(x$VariableUnitsID[i]),
          as.character(x$SampleMedium[i]),
          as.character(x$ValueType[i]),
          as.character(x$IsRegular[i]),
          as.numeric(x$TimeSupport[i]),
          as.numeric(x$TimeUnitsID[i]),
          as.character(x$DataType[i]),
          as.character(x$GeneralCategory[i]),
          as.numeric(x$NoDataValue[i]),
          as.character(status.code)
        )
      added.variables[nrow(added.variables)+1,] <- new.variable
    }
  }
  return (added.variables)
}
