#' GetSiteInfo
#'
#' This function gets the table variables measured at a specific site from the WaterML web service
#'
#' @import XML
#' @param server The URL of the web service ending with .asmx or .wsdl,
#'  for example: http://worldwater.byu.edu/app/index.php/rushvalley/services/cuahsi_1_1.asmx?WSDL
#' @param siteCode The full site code, for example: default:Ru5BMMA. To get a list of
#' available site codes, see GetSites() function and use the FullSiteCode field.
#' @return a data.frame of data values with the following columns:
#' \tabular{lll}{
#' Name \tab Type \tab Description \cr
#' SiteID \tab character \tab The site ID in the original database\cr
#' SiteName \tab character \tab The name of the site \cr
#' SiteCode \tab character \tab A short unique code of the site\cr
#' FullSiteCode \tab character \tab The complete unique code of the site \cr
#'               in the format NETWORK:CODE, for example SNOTEL:879.
#'               \cr
#' Latitude \tab  numeric \tab The WGS84 latitude in decimal degrees \cr
#' Longitude \tab numeric \tab The WGS84 longitude in decimal degrees \cr
#' Elevation \tab numeric \tab The elevation of the site above sea level in meters \cr
#' State \tab character \tab Only for sites in the USA: the state of the site \cr
#' County \tab character \tab Only for sites in the USA: The county of the site \cr
#' Comments \tab character \tab Additional comments about the sites \cr
#'                          (note: this field is often empty)
#'                          \cr
#' VariableCode \tab character \tab Short code of the variable \cr
#' FullVariableCode \tab character \tab The full variable code, for example: SNOTEL:SNWD. \cr
#'                   Use this value as the variableCode parameter in GetValues().
#'                   \cr
#' VariableName \tab character \tab The name of the variable \cr
#' ValueType \tab character \tab the type of observation: \cr
#'            Field Observation or Derived Value
#'            \cr
#' DataType \tab character \tab the aggregate data type: \cr
#'           Average, Continuous, Sporadic..
#'           \cr
#' GeneralCategory \tab character \tab the general category of the measurements: \cr
#'                  Climate, Water Quality..
#'                  \cr
#' SampleMedium \tab character \tab the sample medium: \cr
#'               for example water, atmosphere, soil..
#'               \cr
#' UnitName \tab character \tab The name of the measurement units \cr
#' UnitType \tab character \tab the type of the measurement units \cr
#' UnitAbbreviation \tab character \tab The abbreviation of the measurement units \cr
#'                   (m, cm, in..)
#'                   \cr
#' NoDataValue \tab numeric \tab The value that indicates missing data \cr
#' IsRegular \tab boolean \tab TRUE if the measurements are regular, FALSE otherwise \cr
#' TimeUnitName \tab character \tab The name of the time units \cr
#' TimeUnitAbbreviation \tab character \tab The time units abbreviation \cr
#' TimeSupport \tab character \tab The length of the time period over which \cr
#'              one measurement is taken
#'              \cr
#' Speciation \tab character \tab The chemical sample speciation \cr
#'             (as nitrogen, as phosphorus..)
#'             \cr
#' methodID \tab character \tab The ID of the sensor or measurement method \cr
#' methodCode \tab character \tab The code of the sensor or measurement method. \cr
#'             Usually the same as methodID.
#'             \cr
#' methodDescription \tab character \tab The description of the sensor or \cr
#'                    of the data collection instrumentation / measurement method.
#'                    \cr
#' methodLink \tab character \tab The hyperlink of the website \cr
#'             of the sensor or measurement method.
#'             \cr
#' sourceID \tab character \tab The ID of the data source or author \cr
#' organization \tab character \tab The name of the organization collecting the data \cr
#' sourceDescription \tab character \tab The description of organization collecting the data \cr
#' citation \tab character \tab Instruction how to cite the data \cr
#' qualityControlLevelID \tab character \tab The ID of the quality control level. \cr
#'                        Usually 0 means raw data and 1 means quality controlled data.
#'                        \cr
#' qualityControlLevelCode: \tab character \tab The code of the quality control level.
#'                           Usually same as qualityControlLevelID.
#'                           \cr
#' qualityControlLevelDefinition: \tab character \tab The quality control level definition. \cr
#' valueCount: \tab character \tab The number of observations in this time series \cr
#' beginDateTime: \tab POSIXct \tab The local date and time of the first available \cr
#'                              observation in this time series.
#'                              \cr
#' endDateTime: \tab POSIXct \tab The local date and time of the last available \cr
#'                            observation in this time series.
#'                            \cr
#' beginDateTimeUTC: \tab POSIXct \tab The UTC date and time of the last available
#'                                 observation in this time series.
#'                                 \cr
#' endDateTimeUTC: \tab POSIXct \tab The UTC date and time of the last available
#'                               observation in this time series.
#'                               \cr
#' }
#' The output data.frame also has attributes with information about the status:
#' download.time, parse.time, download.status, parse.status
#' These attributes can be used for troubleshooting WaterOneFlow/WaterML server errors.
#' If parse status is "NO_SERIES_FOUND", then this site doesn't have any available data.
#' @keywords waterml
#' @export
#' @examples
#' server <- "http://worldwater.byu.edu/app/index.php/rushvalley/services/cuahsi_1_1.asmx"
#' siteInfo <- GetSiteInfo(server, siteCode="default:Ru5BMMA")

GetSiteInfo <- function(server, siteCode) {

  # declare the default download timeout in seconds
  max_timeout = 360

  # declare empty return data frame
  df <- data.frame()

  # trim any leading and trailing whitespaces in server
  server <- gsub("^\\s+|\\s+$", "", server)

  # if server ends with ?WSDL or ?wsdl, we assume that service is SOAP
  # otherwise, assume that service is REST
  SOAP <- TRUE

  # if server ends with .asmx, we also assume that the service is SOAP and we add ?WSDL
  m1 <- regexpr("asmx$", server)
  if (m1 > 1) {
    server <- paste(server, "WSDL", sep="?")
  }


  # if server ends with ?WSDL or ?wsdl, we assume that service is SOAP
  # otherwise, assume that service is REST
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
    methodName <- "GetSiteInfoObject"

    SOAPAction <- paste(namespace, methodName, sep="")
    headers <- c("Content-Type" = "text/xml","SOAPAction" = SOAPAction)
    envelope <- MakeSOAPEnvelope(namespace, methodName, c(site=siteCode))

    print(paste("downloading SiteInfo from:", url))

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
    #if the service is REST
    print(paste("downloading SiteInfo from:", server))
    downloaded <- FALSE
    download.time <- system.time(
      err <- tryCatch({
        response <- GET(server, timeout(max_timeout))
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
    version <- "1.1"

    print(paste("download time:", download.time["elapsed"], "seconds, status:", status.code))
  }
  attr(df, "download.time") <- download.time["elapsed"]
  attr(df, "download.status") <- status.code

  ######################################################
  # Parsing the WaterML XML Data                       #
  ######################################################
  begin.parse.time <- Sys.time()

  doc <- xmlParse(response)

  # specify the namespace information
  ns <- WaterOneFlowNamespace(version)

  #try to find faultstring to look for an error
  fault <- xpathSApply(doc, "//soap:Fault", xmlValue, namespaces=ns)
  if (length(fault) > 0) {
    print(paste("SERVER ERROR in GetSiteInfo ", as.character(fault), sep=":"))
    attr(df, "download.time") <- download.time["elapsed"]
    attr(df, "download.status") <- as.character(fault)
    attr(df, "parse.time") <- NA
    attr(df, "parse.status") <- "SERVER FAULT"
    return(df)
  }

  SiteName = xpathSApply(doc, "//sr:siteName", xmlValue, namespaces=ns)
  N <- length(SiteName)
  SiteCode = xpathSApply(doc, "//sr:siteCode", xmlValue, namespaces=ns)
  Network = xpathSApply(doc, "//sr:siteCode", xmlGetAttr, name="network", namespaces=ns)

  SiteID <- xpathSApply(doc, "//sr:siteCode", xmlGetAttr, name="siteID", namespaces=ns)
  SiteID <- unlist(SiteID)
  numSiteIDs <- length(SiteID)
  if (numSiteIDs != N) {
    SiteID <- SiteCode
  }

  Latitude <- xpathSApply(doc, "//sr:latitude", xmlValue, namespaces=ns)
  Longitude = xpathSApply(doc, "//sr:longitude", xmlValue, namespaces=ns)

  Elevation <- xpathSApply(doc, "//sr:elevation_m", xmlValue, namespaces=ns)
  numElevations <- length(Elevation)
  if (numElevations != N) {
    Elevation <- NA
  }

  # State, County, Comments: different tags for WaterML 1.0 and 1.1
  if (version=="1.1"){
    State = xpathSApply(doc, "//sr:siteProperty[@name='State']", xmlValue, namespaces=ns)
    County = xpathSApply(doc, "//sr:siteProperty[@name='County']", xmlValue, namespaces=ns)
    Comments = xpathSApply(doc, "//sr:siteProperty[@name='Site Comments']", xmlValue, namespaces=ns)
  } else {
    State = xpathSApply(doc, "//sr:note[@title='State']", xmlValue, namespaces=ns)
    County = xpathSApply(doc, "//sr:note[@title='County']", xmlValue, namespaces=ns)
    Comments = NA
  }
  # Check for empty values of state, county, comments
  numStates <- length(State)
  if (numStates != N) {
    State <- NA
  }
  numCounties <- length(County)
  if (numCounties != N) {
    County <- NA
  }
  numComments <- length(Comments)
  if (numComments != N) {
    Comments <- NA
  }

  VariableCode <- xpathSApply(doc, "//sr:variableCode", xmlValue, namespaces=ns)
  N <- length(VariableCode)

  # Check for 'No Series Found' case
  if (N==0) {
    print(paste("NOTE: 0 time series found for site:", siteCode))

    end.parse.time <- Sys.time()
    parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))

    attr(df, "download.time") <- download.time["elapsed"]
    attr(df, "download.status") <- "success"
    attr(df, "parse.time") <- parse.time
    attr(df, "parse.status") <- "NO_SERIES_FOUND"
    return(df)
  }

  VariableName <- xpathSApply(doc, "//sr:variableName", xmlValue, namespaces=ns)

  VariableID <- unlist(xpathSApply(doc, "//sr:variableCode", xmlGetAttr, name="variableID", namespaces=ns))
  if (length(VariableID) == 0) { VariableID <- VariableCode }

  Vocabulary <- unlist(xpathSApply(doc, "//sr:variableCode", xmlGetAttr, name="vocabulary", namespaces=ns))

  #######################################################################################
  # START of SPECIAL CASE: process variable: use special case for WaterML 1.0           #
  #for WaterML 1.0 we must use a loop, because elements with unknown values are missing #
  #######################################################################################
  if (version == "1.0") {
    allVariables <- getNodeSet(doc, "//sr:series/sr:variable", namespaces=ns)
    i <- 1
    if (length(allVariables) < N) {
      print("Bad XML format: not enough details about the variables")

      end.parse.time <- Sys.time()
      parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))

      attr(df, "download.time") <- download.time["elapsed"]
      attr(df, "download.status") <- "success"
      attr(df, "parse.time") <- parse.time
      attr(df, "parse.status") <- "BAD_XML_FORMAT"
      return(df)
    }
    #allocate vectors for variables
    ValueType=rep("",N)
    DataType=rep("",N)
    GeneralCategory=rep("",N)
    SampleMedium=rep("",N)
    UnitName=rep("",N)
    UnitType=rep("",N)
    UnitAbbreviation=rep("",N)
    NoDataValue=rep(NA,N)
    IsRegular=rep("",N)
    TimeUnitName=rep("",N)
    TimeUnitAbbreviation=rep("",N)
    TimeSupport=rep("",N)
    Speciation=rep("",N)

    for (i in 1:N) {
      varObj <- allVariables[[i]]
      v <- unlist(xmlToList(varObj))
      ValueType[i] <- v["valueType"]
      DataType[i] <- v["dataType"]
      if (!is.null(v["generalCategory"])) {
        GeneralCategory[i] <- v["generalCategory"]
      }
      if (!is.null(v["sampleMedium"])) {
        SampleMedium[i] <- v["sampleMedium"]
      }
      UnitName[i] <- v["units.text"]
      if (!is.null(v["units..attrs.unitsType"])) {
        UnitType[i] <- v["units..attrs.unitsType"]
      }
      UnitAbbreviation[i] <- v["units..attrs.unitsAbbreviation"]
      IsRegular[i] <- ifelse(is.na(v["timeSupport..attrs.isRegular"]), v["timeSupport.isRegular"],
                                v["timeSupport..attrs.isRegular"])
      TimeUnitName[i] <- v["timeSupport.unit.UnitDescription"]
      if (is.na(TimeUnitName[i])) {
        TimeUnitName[i] <- v["timeSupport.unit.UnitName"]
      }
      TimeUnitAbbreviation[i] <- v["timeSupport.unit.UnitAbbreviation"]
      TimeSupport[i] <- v["timeSupport.timeInterval"]

      if (!is.null(v["NoDataValue"])) {
        NoDataValue[i] <- as.numeric(v["NoDataValue"])
      }
    }

    MethodID <- unlist(xpathSApply(doc, "//sr:Method", xmlGetAttr, name="methodID", namespaces=ns))
    MethodCode <- xpathSApply(doc, "//sr:MethodCode", xmlValue, namespaces=ns)

    MethodDescription <- xpathSApply(doc, "//sr:MethodDescription", xmlValue, namespaces=ns)
    if (length(MethodDescription) < N) { MethodDescription <- NA }

    MethodLink <- xpathSApply(doc, "//sr:MethodLink", xmlValue, namespaces=ns)
    if (length(MethodLink) < N) { MethodLink <- NA }

    if (length(MethodID) < N & length(MethodCode) == N) {
      MethodID <- MethodCode
    }
    if (length(MethodCode) < N & length(MethodID) == N) {
      MethodCode <- MethodID
    }
    if (length(MethodID) < N) { MethodID <- NA }
    if (length(MethodCode) < N) { MethodCode <- NA }

    SourceID <- unlist(xpathSApply(doc, "//sr:Source", xmlGetAttr, name="sourceID", namespaces=ns))
    if (length(SourceID) < N) { SourceID <- NA }

    Organization <- xpathSApply(doc, "//sr:Organization", xmlValue, namespaces=ns)
    if (length(Organization) < N) { Organization <- NA }

    SourceDescription <- xpathSApply(doc, "//sr:SourceDescription", xmlValue, namespaces=ns)
    if (length(SourceDescription) < N) { SourceDescription <- NA }

    Citation <- xpathSApply(doc, "//sr:Citation", xmlValue, namespaces=ns)
    if (length(Citation) < N) { Citation <- NA }

    QualityControlLevelID <- unlist(xpathSApply(doc, "//sr:QualityControlLevel", xmlGetAttr,
                                      name="qualityControlLevelID", namespaces=ns))

    QualityControlLevelCode <- xpathSApply(doc, "//sr:qualityControlLevelCode", xmlValue, namespaces=ns)

    if (length(QualityControlLevelID) < N & length(QualityControlLevelCode == N)) {
      QualityControlLevelID <- QualityControlLevelCode
    }
    if (length(QualityControlLevelCode) < N & length(QualityControlLevelID == N)) {
      QualityControlLevelCode <- QualityControlLevelID
    }
    if (length(QualityControlLevelID) < N) { QualityControlLevelID <- NA }
    if (length(QualityControlLevelCode) < N) { QualityControlLevelCode <- NA }

    QualityControlLevelDefinition=xpathSApply(doc, "//sr:QualityControlLevel", xmlValue, namespaces=ns)
    if (length(QualityControlLevelDefinition) < N) { QualityControlLevelDefinition <- NA }

    ValueCount <- xpathSApply(doc, "//sr:valueCount", xmlValue, namespaces=ns)

    BeginDateTime <- xpathSApply(doc, "//sr:beginDateTime", xmlValue, namespaces=ns)
    EndDateTime <- xpathSApply(doc, "//sr:endDateTime", xmlValue, namespaces=ns)

    BeginDateTimeUTC <- xpathSApply(doc, "//sr:beginDateTimeUTC", xmlValue, namespaces=ns)
    if (length(BeginDateTimeUTC) == 0) { BeginDateTimeUTC <- BeginDateTime }

    EndDateTimeUTC <- xpathSApply(doc, "//sr:endDateTimeUTC", xmlValue, namespaces=ns)
    if (length(EndDateTimeUTC) == 0) { EndDateTimeUTC <- EndDateTime }
  #################################################################################################
  # END of SPECIAL CASE of WaterML 1.0                                                            #
  #################################################################################################
  } else {

    ValueType <- xpathSApply(doc, "//sr:valueType", xmlValue, namespaces=ns)
    DataType <- xpathSApply(doc, "//sr:dataType", xmlValue, namespaces=ns)
    GeneralCategory <- xpathSApply(doc, "//sr:generalCategory", xmlValue, namespaces=ns)
    SampleMedium <- xpathSApply(doc, "//sr:sampleMedium", xmlValue, namespaces=ns)

    UnitName <- xpathSApply(doc, "//sr:units/sr:unitName", xmlValue, namespaces=ns)
    UnitType <- xpathSApply(doc, "//sr:units/sr:unitType", xmlValue, namespaces=ns)
    UnitAbbreviation <- xpathSApply(doc,
      "//sr:variable/sr:units/*[self::sr:unitsAbbreviation or self::sr:unitAbbreviation]",
      xmlValue, namespaces=ns)

    #if UnitName is not found, then we look for /variable/unit instead
    if (length(UnitName) == 0) {
      UnitName <- xpathSApply(doc, "//sr:variable/sr:unit/sr:unitName", xmlValue, namespaces=ns)
      UnitType <- xpathSApply(doc, "//sr:variable/sr:unit/sr:unitType", xmlValue, namespaces=ns)
      UnitAbbreviation <- xpathSApply(doc,
        "//sr:variable/sr:unit/*[self::sr:unitsAbbreviation or self::sr:unitAbbreviation]",
        xmlValue, namespaces=ns)
    }

    NoDataValue <- xpathSApply(doc, "//sr:noDataValue", xmlValue, namespaces=ns)

    IsRegular <- unlist(xpathSApply(doc, "//sr:timeScale", xmlGetAttr, name="isRegular", namespaces=ns))
    if (length(IsRegular) < N) {
      IsRegular <- (DataType != "Sporadic")
    }

    TimeUnitName <- xpathSApply(doc, "//sr:timeScale/sr:unit/sr:unitName", xmlValue, namespaces=ns)
    TimeUnitName <- unlist(TimeUnitName)
    if (length(TimeUnitName) < N) { TimeUnitName <- NA }

    TimeUnitAbbreviation <- xpathSApply(doc,
      "//sr:timeScale/sr:unit/*[self::sr:unitsAbbreviation or self::sr:unitAbbreviation]", xmlValue, namespaces=ns)
    TimeUnitAbbreviation <- unlist(TimeUnitAbbreviation)
    if (length(TimeUnitAbbreviation) < N) { TimeUnitAbbreviation <- NA }

    TimeSupport <- xpathSApply(doc, "//sr:timeSupport", xmlValue, namespaces=ns)
    Speciation <- xpathSApply(doc, "//sr:variable/sr:speciation", xmlValue, namespaces=ns)

    BeginDateTime <- xpathSApply(doc, "//sr:beginDateTime", xmlValue, namespaces=ns)
    EndDateTime <- xpathSApply(doc, "//sr:endDateTime", xmlValue, namespaces=ns)

    BeginDateTimeUTC <- xpathSApply(doc, "//sr:beginDateTimeUTC", xmlValue, namespaces=ns)
    if (length(BeginDateTimeUTC) == 0) { BeginDateTimeUTC <- BeginDateTime }

    EndDateTimeUTC <- xpathSApply(doc, "//sr:endDateTimeUTC", xmlValue, namespaces=ns)
    if (length(EndDateTimeUTC) == 0) { EndDateTimeUTC <- EndDateTime }

    ValueCount <- xpathSApply(doc, "//sr:valueCount", xmlValue, namespaces=ns)

    MethodID <- unlist(xpathSApply(doc, "//sr:method", xmlGetAttr, name="methodID", namespaces=ns))

    MethodCode <- xpathSApply(doc, "//sr:methodCode", xmlValue, namespaces=ns)

    MethodDescription <- xpathSApply(doc, "//sr:methodDescription", xmlValue, namespaces=ns)
    if (length(MethodDescription) < N) { MethodDescription <- NA }

    MethodLink <- xpathSApply(doc, "//sr:methodLink", xmlValue, namespaces=ns)
    if (length(MethodLink) < N) { MethodLink <- NA }

    if (length(MethodID) < N & length(MethodCode) == N) {
      MethodID <- MethodCode
    }
    if (length(MethodCode) < N & length(MethodID) == N) {
      MethodCode <- MethodID
    }
    if (length(MethodID) < N) { MethodID <- NA }
    if (length(MethodCode) < N) { MethodCode <- NA }

    SourceID <- unlist(xpathSApply(doc, "//sr:source", xmlGetAttr, name="sourceID", namespaces=ns))
    if (length(SourceID) < N) { SourceID <- NA }

    Organization <- xpathSApply(doc, "//sr:organization", xmlValue, namespaces=ns)
    if (length(Organization) < N) { Organization <- NA }

    SourceDescription <- xpathSApply(doc, "//sr:sourceDescription", xmlValue, namespaces=ns)
    if (length(SourceDescription) < N) { SourceDescription <- NA }

    Citation <- xpathSApply(doc, "//sr:citation", xmlValue, namespaces=ns)
    if (length(Citation) < N) { Citation <- NA }

    QualityControlLevelID=unlist(xpathSApply(doc, "//sr:qualityControlLevel", xmlGetAttr,
                                      name="qualityControlLevelID", namespaces=ns))
    QualityControlLevelCode=xpathSApply(doc, "//sr:qualityControlLevelCode", xmlValue, namespaces=ns)

    if (length(QualityControlLevelID) < N & length(QualityControlLevelCode == N)) {
      QualityControlLevelID <- QualityControlLevelCode
    }
    if (length(QualityControlLevelCode) < N & length(QualityControlLevelID == N)) {
      QualityControlLevelCode <- QualityControlLevelID
    }
    if (length(QualityControlLevelID) < N) { QualityControlLevelID <- NA }
    if (length(QualityControlLevelCode) < N) { QualityControlLevelCode <- NA }

    QualityControlLevelDefinition=xpathSApply(doc, "//sr:definition", xmlValue, namespaces=ns)
    if (length(QualityControlLevelDefinition) < N) { QualityControlLevelDefinition <- NA }
  }

  #define the columns for the output data frame
  df <- data.frame(SiteName=rep(SiteName, N),
                   SiteID=rep(SiteID, N),
                   SiteCode=rep(SiteCode, N),
                   FullSiteCode = rep(paste(Network, SiteCode, sep=":"), N),
                   Latitude=rep(as.numeric(Latitude), N),
                   Longitude=rep(as.numeric(Longitude), N),
                   Elevation=rep(as.numeric(Elevation), N),
                   State=rep(State, N),
                   County=rep(County, N),
                   Comments=rep(Comments, N),
                   VariableCode=VariableCode,
                   FullVariableCode=paste(Vocabulary, VariableCode, sep=":"),
                   VariableName=VariableName,
                   ValueType=ValueType,
                   DataType=DataType,
                   GeneralCategory=GeneralCategory,
                   SampleMedium=SampleMedium,
                   UnitName=UnitName,
                   UnitType=UnitType,
                   UnitAbbreviation=UnitAbbreviation,
                   NoDataValue=as.numeric(NoDataValue),
                   IsRegular=IsRegular,
                   TimeUnitName=TimeUnitName,
                   TimeUnitAbbreviation=TimeUnitAbbreviation,
                   TimeSupport=TimeSupport,
                   Speciation=Speciation,
                   methodID=MethodID,
                   methodCode=MethodCode,
                   methodDescription=MethodDescription,
                   methodLink=MethodLink,
                   sourceID=SourceID,
                   organization=Organization,
                   sourceDescription=SourceDescription,
                   citation=Citation,
                   qualityControlLevelID=QualityControlLevelID,
                   qualityControlLevelCode=QualityControlLevelCode,
                   qualityControlLevelDefinition=QualityControlLevelDefinition,
                   valueCount=ValueCount,
                   beginDateTime=as.POSIXct(strptime(BeginDateTime, "%Y-%m-%dT%H:%M:%S")),
                   endDateTime=as.POSIXct(strptime(EndDateTime, "%Y-%m-%dT%H:%M:%S")),
                   beginDateTimeUTC=as.POSIXct(strptime(BeginDateTimeUTC, "%Y-%m-%dT%H:%M:%S")),
                   endDateTimeUTC=as.POSIXct(strptime(EndDateTimeUTC, "%Y-%m-%dT%H:%M:%S")),
                   stringsAsFactors=FALSE)

  if (nrow(df) == 0) {
    print(paste("NOTE: 0 time series found for site:", siteCode))

    end.parse.time <- Sys.time()
    parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))

    attr(df, "download.time") <- download.time["elapsed"]
    attr(df, "download.status") <- "success"
    attr(df, "parse.time") <- parse.time
    attr(df, "parse.status") <- "NO_SERIES_FOUND"
    return(df)
  }

  end.parse.time <- Sys.time()
  parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
  attr(df, "download.time") <- download.time["elapsed"]
  attr(df, "download.status") <- "success"
  attr(df, "parse.time") <- parse.time
  attr(df, "parse.status") <- "OK"
  return(df)
}
