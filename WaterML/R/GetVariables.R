#' GetVariables
#'
#' This function gets the table of variables from the WaterML web service
#'
#' @import XML
#' @param server The URL of the web service ending with ?WSDL,
#'  for example: http://worldwater.byu.edu/app/index.php/rushvalley/services/cuahsi_1_1.asmx?WSDL
#' @return a data.frame of variables with the following columns:
#' \tabular{lll}{
#' VariableCode \tab character \tab Short code of the variable \cr
#' FullVariableCode \tab character \tab The full variable code, for example: SNOTEL:897. Use this value
#'                   as the variableCode parameter in GetValues() function.
#'                   \cr
#' VariableName \tab character \tab The name of the variable \cr
#' ValueType \tab character \tab the type of observation: Field Observation or Derived Value \cr
#' DataType \tab character \tab the aggregate data type: Average, Continuous, Sporadic.. \cr
#' GeneralCategory \tab character \tab the general category of the measurements: Climate, Water Quality.. \cr
#' SampleMedium \tab character \tab the sample medium, for example water, atmosphere, soil.. \cr
#' UnitName \tab character \tab The name of the measurement units \cr
#' UnitType \tab character \tab the type of the measurement units \cr
#' UnitAbbreviation \tab character \tab The abbreviation of the measurement units (m, cm, in..) \cr
#' NoDataValue \tab numeric \tab The value that indicates missing data \cr
#' IsRegular \tab boolean \tab TRUE if the measurements are regular, FALSE otherwise \cr
#' TimeUnitName \tab character \tab The name of the time units \cr
#' TimeUnitAbbreviation \tab character \tab The time units abbreviation \cr
#' TimeSupport \tab character \tab The length of the time period over which one measurement is taken \cr
#' Speciation \tab character \tab The chemical sample speciation (as nitrogen, as phosphorus..) \cr
#' }
#' @keywords WaterML
#' @export
#' @examples
#' GetVariables("http://worldwater.byu.edu/app/index.php/rushvalley/services/cuahsi_1_1.asmx?WSDL")

GetVariables <- function(server) {

  # declare the default download timeout in seconds
  max_timeout = 360

  # declare empty return data frame
  df <- data.frame()

  # trim any leading and trailing whitespaces in server
  server <- gsub("^\\s+|\\s+$", "", server)

  # if server ends with .asmx, we also assume that the service is SOAP and we add ?WSDL
  m1 <- regexpr("asmx$", server)
  if (m1 > 1) {
    server <- paste(server, "WSDL", sep="?")
  }

  # if server ends with ?WSDL or ?wsdl, we assume that service is SOAP
  # otherwise, assume that service is REST
  SOAP <- TRUE
  m <- regexpr("?WSDL|wsdl", server)
  if (m > 1) {
    url <- substr(server, 0, m - 2)
    SOAP <- TRUE
  } else {
    SOAP <- FALSE
  }

  #if the service is SOAP:
  if (SOAP) {
    versionInfo <- WaterOneFlowVersion(server)
    namespace <- versionInfo$Namespace
    version <- versionInfo$Version
    methodName <- "GetVariableInfoObject"

    SOAPAction <- paste(namespace, methodName, sep="")
    envelope <- MakeSOAPEnvelope(namespace, methodName, c(variable=""))
    headers <- c("Content-Type" = "text/xml", "SOAPAction" = SOAPAction)

    print(paste("GetVariables from", url))

    downloaded <- FALSE
    download.time <- system.time(
      err <- tryCatch({
        response <- POST(url, body = envelope, add_headers(headers),
                         timeout(max_timeout))
        status <- http_status(response)$message
        downloaded <- TRUE
      },error = function(e) {
        print(conditionMessage(e))
      }
      )
    )
    if (!downloaded) {
      attr(df, "download.time") <- download.time["elapsed"]
      attr(df, "download.status") <- err
      attr(df, "parse.time") <- NA
      attr(df, "parse.status") <- NA
      return(df)
    }

    status.code <- http_status(response)$category
    print(paste("download time:", download.time["elapsed"], "seconds, status:", status.code))

  } else {
    #GetVariables using REST
    print(paste("downloading variables from:", server, "..."))

    downloaded <- FALSE
    download.time <- system.time(
      err <- tryCatch({
        response <- GET(server)
        status <- http_status(response)$message
        downloaded <- TRUE
      },error=function(e){
        print(conditionMessage(e))
      })
    )

    if (!downloaded) {
      attr(df, "download.time") <- download.time["elapsed"]
      attr(df, "download.status") <- err
      attr(df, "parse.time") <- NA
      attr(df, "parse.status") <- NA
      return(df)
    }

    status.code <- http_status(response)$category
    print(paste("download time:", download.time["elapsed"], "seconds, status:", status.code))
  }

  # check for server error category
  if (http_status(response)$category == "server error") {
    attr(df, "download.time") <- download.time["elapsed"]
    attr(df, "download.status") <- http_status(response)$message
    attr(df, "parse.time") <- NA
    attr(df, "parse.status") <- NA
    return(df)
  }


  attr(df, "download.time") <- download.time["elapsed"]
  attr(df, "download.status") <- status.code

  ######################################################
  # Parsing the WaterML XML Data                       #
  ######################################################

  begin.parse.time <- Sys.time()

  print("reading variables WaterML data...")
  doc <- NULL
  err <- tryCatch({
    doc <- xmlParse(response)
  }, warning = function(w) {
    print("Error reading WaterML: Bad XML format.")
    attr(df, "parse.status") <- "Bad XML format"
    attr(df, "parse.time") <- 0
    return(df)
  }, error = function(e) {
    print("Error reading WaterML: Bad XML format.")
    attr(df, "parse.status") <- "Bad XML format"
    attr(df, "parse.time") <- 0
    return(df)
  }
  )
  if (is.null(doc)) {
    print("Error reading WaterML: Bad XML format.")
    attr(df, "parse.status") <- "Bad XML format"
    attr(df, "parse.time") <- 0
    return(df)
  }

  # specify the namespace information
  ns <- WaterOneFlowNamespace(version)

  vars <- getNodeSet(doc, "//sr:variable", namespaces=ns)

  N <- xmlSize(vars)
  #define the columns
  df <- data.frame(
    VariableCode=rep("",N), FullVariableCode=rep("",N),
    VariableName=rep("",N), ValueType=rep("",N),
    DataType=rep("",N), GeneralCategory=rep("",N),SampleMedium=rep("",N),
    UnitName=rep(NA,N), UnitType=rep(NA,N), UnitAbbreviation=rep(NA,N),
    NoDataValue=rep(NA,N), IsRegular=rep("",N),
    TimeUnitName=rep("",N), TimeUnitAbbreviation=rep("",N),
    TimeSupport=rep("",N), Speciation=rep("",N),
    stringsAsFactors=FALSE)

  for(i in 1:N) {
    varObj <- vars[[i]]
    v <- unlist(xmlToList(varObj))
    varcode <- v["variableCode.text"]
    df$VariableCode[i] <- varcode
    df$FullVariableCode[i] <- paste(v["variableCode..attrs.vocabulary"], varcode, sep=":")
    df$VariableName[i] <- v["variableName"]
    df$ValueType[i] <- v["valueType"]
    df$DataType[i] <- v["dataType"]
    df$GeneralCategory[i] <- v["generalCategory"]
    df$SampleMedium[i] <- v["sampleMedium"]
    if (version == "1.1") {
      df$UnitName[i] <- ifelse(is.na(v["unit.unitName"]), v["units.unitName"], v["unit.unitName"])
      df$UnitType[i] <- v["unit.unitType"]
      df$UnitAbbreviation[i] <- v["unit.unitAbbreviation"]
      df$IsRegular[i] <- ifelse(is.na(v["timeScale..attrs.isRegular"]), v["timeScale.isRegular"],
                             v["timeScale..attrs.isRegular"])

      df$TimeUnitName[i] <- v["timeScale.unit.unitName"]
      df$TimeUnitAbbreviation[i] <- v["timeScale.unit.unitAbbreviation"]
      df$TimeSupport[i] <- v["timeScale.timeSupport"]
      df$NoDataValue[i] <- as.numeric(v["noDataValue"])
    } else {
      df$UnitName[i] <- v["units.text"]
      df$UnitType[i] <- v["units..attrs.unitsType"]
      df$UnitAbbreviation[i] <- v["units..attrs.unitsAbbreviation"]
      df$IsRegular[i] <- ifelse(is.na(v["timeSupport..attrs.isRegular"]), v["timeSupport.isRegular"],
                                   v["timeSupport..attrs.isRegular"])
      df$TimeUnitName[i] <- v["timeSupport.unit.UnitDescription"]
      df$TimeUnitAbbreviation[i] <- v["timeSupport.unit.UnitAbbreviation"]
      df$TimeSupport[i] <- v["timeSupport.timeInterval"]
      df$NoDataValue[i] <- as.numeric(v["NoDataValue"])
    }
    df$Speciation[i] <- v["speciation"]
  }

  end.parse.time <- Sys.time()
  parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
  attr(df, "download.time") <- download.time["elapsed"]
  attr(df, "download.status") <- "success"
  attr(df, "parse.time") <- parse.time
  attr(df, "parse.status") <- "OK"
  return(df)
}
