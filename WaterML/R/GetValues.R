#' GetValues
#'
#' This function gets the time series data values from the WaterML web service
#'
#' @import stats
#' @import XML
#' @import httr
#' @param server The URL of the web service,
#'  for example: http://worldwater.byu.edu/interactive/rushvalley/services/index.php/cuahsi_1_1.asmx?WSDL.
#'  This can be also a custom REST URL or the file name of the WaterML file.
#' @param siteCode The site code. To get a list of available site codes, see GetSites() function
#'  and use the FullSiteCode field.
#' @param variableCode The variable code. To get a list of possible variable codes, see GetVariables()
#'  function and use the FullVariableCode field
#' @param startDate (optional) The start date in "yyyy-mm-dd" format
#' @param endDate (optional) The end date in "yyyy-mm-dd" format
#' @param methodID (optional) The ID of the observation method. To get a list of possible method IDs, see
#' methodID column in the output of GetSiteInfo(). If methodID is not specified, then the observations
#' in the output data.frame won't be filtered by method.
#' @param sourceID (optional) The ID of the source. To get a list of possible source IDs, see
#' sourceID column in the output of GetSiteInfo(). If sourceID is not specified, then the observations
#' in the output data.frame won't be filtered by source.
#' @param qcID (optional) The ID of the quality control level. Typically 0 is used for raw data and 1 is
#' used for quality controlled data. To get a list of possible quality control level IDs, see
#' QualityControlLevelID column in the output of GetSiteInfo(). If qcID is not specified, then the
#' observations in the output data.frame won't be filtered by quality control level.
#' @param daily (optional) If you set daily="max", daily="min" or daily="mean", then the
#' data values are aggreagted to daily time step.
#' @return a data.frame of data values with the following columns:
#' \itemize{
#' \item time: The local date/time of the observation. The data type is POSIXct. POSIXct is
#'             a data type in R for storing time.
#' \item DataValue: The observed data value
#' \item UTCOffset: The difference between local time and UTC time in hours
#' \item CensorCode: The code for censored observations. Possible values are nc (not censored),
#'             gt (greater than), lt (less than),
#'             nd (non-detect), pnq (present but not quantified)
#' \item DateTimeUTC: The UTC time of the observation. The data type is POSIXct.
#'             POSIXct is a special data type in R for storing time.
#' \item MethodCode: The code of the method or instrument used for the observation
#' \item SourceCode: The code of the data source
#' \item QualityControlLevelCode: The code of the quality control level. Possible values are
#'             -9999 (Unknown), 0 (Raw data), 1 (Quality controlled data),
#'             2 (Derived products), 3 (Interpreted products), 4 (Knowledge products)
#' }
#' The output data.frame also has attributes with information about the status:
#' download.time, parse.time, download.status, parse.status
#' These attributes can be used for troubleshooting WaterOneFlow/WaterML server errors.
#' If parse status is "NO_VALUES_FOUND",
#' then this time series doesn't have any available data for the selected time period.
#' @keywords waterml
#' @export
#' @examples
#' #example 1: Get Values from a known site and variable from RushValley server
#' v1 <- GetValues("http://worldwater.byu.edu/app/index.php/rushvalley/services/cuahsi_1_1.asmx?WSDL",
#'            site="Ru5BMMA", variable="SRS_Nr_NDVI", startDate="2014-11-01", endDate="2014-11-02",
#'            daily="max")
#' #example 2: Get values from an external REST URL (in this case the Provo USGS NWIS site id 10163000)
#' url <- "http://waterservices.usgs.gov/nwis/dv/?format=waterml,1.1&sites=10163000&parameterCd=00060"
#' v2 <- GetValues(url)
#' #example 3: Get values from WaterML 2.0 file and show year, month, day
#' url2 <- "http://www.waterml2.org/KiWIS-WML2-Example.wml"
#' waterml2_data <- GetValues(url2)
#' waterml2_data$year <- strftime(waterml2_data$time, "%Y")
#' waterml2_data$month <- strftime(waterml2_data$time, "%M")
#' waterml2_data$day <- strftime(waterml2_data$time, "%d")

GetValues <- function(server, siteCode=NULL, variableCode=NULL, startDate=NULL, endDate=NULL,
                      methodID=NULL, sourceID=NULL, qcID=NULL, daily=NULL) {

  #file or url?
  isFile <- FALSE

  # declare the default download timeout in seconds
  max_timeout = 360

  # declare empty return data frame
  df <- data.frame()

  # trim any leading and trailing whitespaces in server
  server <- gsub("^\\s+|\\s+$", "", server)

  #if server is a file name

  SOAP <- TRUE

  # if server ends with .asmx, we also assume that the service is SOAP and we add ?WSDL
  m1 <- regexpr("asmx$", server)
  if (m1 > 1) {
    server <- paste(server, "WSDL", sep="?")
  }

  #save variableCode for possible future use
  original_variable_code <- NULL
  if (!is.null(variableCode)) {
    original_variable_code <- variableCode
  }

  #check startDate, endDate if it is null
  startDateParam <- ifelse(is.null(startDate), "", strftime(as.POSIXct(startDate), "%Y-%m-%dT%H:%M:%S"))
  endDateParam <- ifelse(is.null(endDate), "", strftime(as.POSIXct(endDate), "%Y-%m-%dT%H:%M:%S"))

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
    methodName <- "GetValuesObject"

    #format the variable with the methodID, sourceID, qcID parameters
    variableCodeParam <- variableCode
    if (!is.null(methodID)) {
      variableCodeParam <- paste(variableCodeParam, ":methodCode=", methodID, sep="")
    }
    if (!is.null(sourceID)) {
      variableCodeParam <- paste(variableCodeParam, ":sourceCode=", sourceID, sep="")
    }
    if (!is.null(qcID)) {
      variableCodeParam <- paste(variableCodeParam, ":qualityControlLevelCode=", qcID, sep="")
    }

    SOAPAction <- paste(namespace, methodName, sep="")
    envelope <- MakeSOAPEnvelope(namespace, methodName, c(location=siteCode,
                                                          variable=variableCodeParam,
                                                          startDate=startDateParam,
                                                          endDate=endDateParam))
    headers <- c("Content-Type" = "text/xml", "SOAPAction" = SOAPAction)

    print(paste("downloading values from:", url, "..."))

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
      attr(df, "download.time") <- as.numeric(download.time["elapsed"])
      attr(df, "download.status") <- err
      attr(df, "parse.time") <- NA
      attr(df, "parse.status") <- NA
      return(df)
    }

    status.code <- http_status(response)$category
    print(paste("download time:", download.time["elapsed"], "seconds, status:", status.code))
    # check for bad status code
    if (tolower(status.code) != "success") {
      status.message <- http_status(response)$message
      attr(df, "download.time") <- as.numeric(download.time["elapsed"])
      attr(df, "download.status") <- status.message
      attr(df, "parse.time") <- NA
      attr(df, "parse.status") <- NA
      return(df)
    }
  } else {
    #REST
    version <- "1.1"

    if (substr(server, 1, 4) == "http")
    {
      print(paste("downloading values from:", server, "..."))

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
        attr(df, "download.time") <- as.numeric(download.time["elapsed"])
        attr(df, "download.status") <- err
        attr(df, "parse.time") <- NA
        attr(df, "parse.status") <- NA
        return(df)
      }

      status.code <- http_status(response)$category
      print(paste("download time:", download.time["elapsed"], "seconds, status:", status.code))
      # check for bad status code
      if (tolower(status.code) != "success") {
        status.message <- http_status(response)$message
        attr(df, "download.time") <- as.numeric(download.time["elapsed"])
        attr(df, "download.status") <- status.message
        attr(df, "parse.time") <- NA
        attr(df, "parse.status") <- status.message
        return(df)
      }
    } else {
      #we are using a local file..
      isFile <- TRUE
    }
  }

  if (!isFile) {
    download.time <- as.numeric(download.time["elapsed"])
    download.status <- status.code
    attr(df, "download.time") <- download.time
    attr(df, "download.status") <- download.status
  } else {
    download.time <- 0
    download.status <- "success"
  }

  ######################################################
  # Parsing the WaterML XML Data                       #
  ######################################################
  begin.parse.time <- Sys.time()

  print("reading data values WaterML ...")
  doc <- NULL
  status.code <- "xml error"
  err <- tryCatch({
    if (isFile) {
      doc <- xmlParseDoc(server)
      status.code <- "success"
    } else {
      doc <- xmlParse(content(response, type="text", encoding="utf-8"))
    }
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
  })
  if (is.null(doc)) {
    print("Error reading WaterML: Bad XML format.")
    attr(df, "parse.status") <- "Bad XML format"
    attr(df, "parse.time") <- 0
    return(df)
  }

  # Check for WaterML version 2.0 (special code - adopted from dataRetrieval package..)
  waterml_version <- WaterMLVersion(doc)
  if (waterml_version == "2.0") {

    ns <- xmlNamespaceDefinitions(doc, simplify = TRUE)
    timeSeries <- xpathApply(doc, "//wml2:Collection", namespaces = ns)

    if(0 == length(timeSeries)){
      df <- data.frame()
      print("NOTE: No data values found in this time series")
      end.parse.time <- Sys.time()
      parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
      attr(df, "parse.status") <- "NO_VALUES_FOUND"
      attr(df, "parse.time") <- parse.time
      return(df)
    }

    for (i in 1:length(timeSeries)){

      chunk <- xmlDoc(timeSeries[[i]])
      chunk <- xmlRoot(chunk)
      chunkNS <- xmlNamespaceDefinitions(chunk, simplify = TRUE)

      xp <- xpathApply(chunk, "//wml2:MeasurementTimeseries/wml2:point/wml2:MeasurementTVP",
                       xpathSApply, ".//*[not(*)]",
                       function(x) setNames(ifelse(nzchar(xmlValue(x)),
                                                   xmlValue(x),
                                                   ifelse("qualifier" == xmlName(x),
                                                          xpathSApply(x,"./@xlink:title",namespaces = ns),"")), #originally I had the "" as xmlAttr(x)
                                            xmlName(x,full=TRUE)),
                       namespaces = chunkNS)

      if(length(xpathApply(doc,
                           "//wml2:MeasurementTimeseries/wml2:point/wml2:MeasurementTVP/wml2:metadata/wml2:TVPMeasurementMetadata",
                           xmlValue, namespaces = ns)) != 0){
        xp <- xp[-1]
      }

      xp2 <- unlist(xp)
      xpTimes <- xp2[seq(1, length(xp2), 2)]
      xpVals <- xp2[seq(2, length(xp2), 2)]

      DF2 <- data.frame(time=xpTimes, value=xpVals, stringsAsFactors = FALSE)

      DF2$time <- substr(gsub(":","",DF2$time),1, 17)
      DF2$time <- ifelse(nchar(DF2$time) > 18,
                         as.POSIXct(DF2$time, format="%Y-%m-%dT%H%M%S%z",tz="UTC"),
                         as.POSIXct(DF2$time, format="%Y-%m-%dT%H%M%S",tz="UTC"))

      DF2$time <- as.POSIXct(DF2$time, origin = "1970-01-01", tz="UTC")
      DF2$value <- as.numeric(DF2$value)

      #########################################
      # Very specific to USGS:
      defaultQualifier <- as.character(xpathApply(chunk, "//wml2:defaultPointMetadata/wml2:DefaultTVPMeasurementMetadata/wml2:qualifier/@xlink:title",namespaces = chunkNS))

      if (length(defaultQualifier) == 0 && (typeof(defaultQualifier) == "character")) {
        defaultQualifier <- "NA"
      }

      if("swe:value" %in% names(DF2)){
        isQual <- as.character(xpathApply(chunk,
                                          "//wml2:MeasurementTimeseries/wml2:point/wml2:MeasurementTVP/wml2:metadata/wml2:TVPMeasurementMetadata/wml2:qualifier/@xlink:title",
                                          namespaces = chunkNS))
        DF2$qualifier <- ifelse(defaultQualifier != isQual,isQual,defaultQualifier)
        DF2$`swe:value` <- NULL
      } else {
        DF2$qualifier <- rep(defaultQualifier,nrow(DF2))
      }

      names(DF2) <- c("time", "DataValue", "Qualifier")
      DF2$UTCOffset <- 0
      DF2$CensorCode <- "nc"
      DF2$DateTimeUTC <- DF2$LocalDateTime
      DF2$MethodCode <- NA
      DF2$SourceCode <- NA
      DF2$QualityControlLevelCode <- NA

      end.parse.time <- Sys.time()
      parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
      attr(DF2, "download.time") <- download.time
      attr(DF2, "download.status") <- download.status
      attr(DF2, "parse.status") <- "OK"
      attr(DF2, "parse.time") <- parse.time
    }
    return (DF2)
  }

  # specify the namespace information
  ns <- WaterOneFlowNamespace(version)

  #try to find faultstring to look for an error
  fault <- xpathSApply(doc, "//soap:Fault", xmlValue, namespaces=ns)
  if (length(fault) > 0) {
    print(paste("SERVER ERROR in GetValues ", as.character(fault), sep=":"))
    end.parse.time <- Sys.time()
    parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
    attr(df, "parse.status") <- as.character(fault)
    attr(df, "parse.time") <- parse.time
    return(df)
  }

  #again check for the status code
  if (tolower(status.code) == "server error") {
    print(paste("SERVER ERROR in GetValues ", http_status(response)$message))
    end.parse.time <- Sys.time()
    parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
    attr(df, "parse.status") <- http_status(response)$message
    attr(df, "parse.time") <- parse.time
    return(df)
  }


  # extract the data columns with XPath
  val = xpathSApply(doc, "//sr:value", xmlValue, namespaces=ns)
  N <- length(val)

  # if N is 0: the document does not have data values or the xml is probably not valid
  if (N == 0) {
    timeSeriesElement <- unlist(xpathSApply(doc, "//sr:timeSeries", xmlValue, namespaces=ns))
    if (is.null(timeSeriesElement)) {
      print("Error reading WaterML: Bad XML format. <timeSeries> tag not found.")
      end.parse.time <- Sys.time()
      parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
      attr(df, "parse.status") <- "Bad XML format. <timeSeries> tag not found."
      attr(df, "parse.time") <- parse.time
      return(df)
    } else {
      #no data values were found
      print("NOTE: No data values found in this time series")

      #special case: methodID, sourceID or qcID is specified. Try again with
      #empty methodID, sourceID, qcID
      if (!is.null(methodID) | !is.null(sourceID) | !is.null(qcID)) {
        print("Trying GetValues again without methodID, sourceID, qcID...")
        daily_param <- daily
        return(GetValues(server, siteCode, original_variable_code, startDate, endDate,
                         methodID=NULL, sourceID=NULL, qcID=NULL, daily=daily_param))
      }

      end.parse.time <- Sys.time()
      parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
      attr(df, "parse.status") <- "NO_VALUES_FOUND"
      attr(df, "parse.time") <- parse.time
      return(df)
    }
  }

  #look for zoneOffset
  time_diff <- NULL
  zoneOffset <- xpathSApply(doc, "//sr:defaultTimeZone", xmlGetAttr, name="zoneOffset", namespaces=ns)
  zoneOffset <- unlist(zoneOffset)
  zoneName <- "GMT"
  if (length(zoneOffset) > 0) {
    offset_split <- strsplit(zoneOffset, ":")
    diff_text <- offset_split[[1]][1]
    time_diff <- as.difftime(as.numeric(diff_text), units="hours")
    if (as.numeric(diff_text) > 0) {
      zoneName <- paste("Etc/GMT+", as.numeric(diff_text), sep="")
    } else {
      zoneName <- paste("Etc/GMT", as.numeric(diff_text), sep="")
    }
  }

  bigData <- 10000
  if (N > bigData) {
    print(paste("found", N,"data values"))
  }

  if (N > bigData) { print("processing censorCode...") }
  censorCode = xpathSApply(doc, "//sr:value", xmlGetAttr, name="censorCode", namespaces=ns)
  censorCode <- unlist(censorCode)
  if (is.null(censorCode)) {
    censorCode <- rep("nc", N)
  }
  if (N > bigData) { print("processing qualifiers...") }
  qualifier <- xpathSApply(doc, "//sr:value", xmlGetAttr, name="qualifiers", namespaces=ns)
  qualifier <- unlist(qualifier)
  if (is.null(qualifier)) {
    qualifier <- rep(NA, N)
  }
  # currently we require that all values have a qualifier attached to it,
  # or none of the values have a qualifier
  if (length(qualifier) < N) {
    qualifier <- rep(NA, N)
  }

  if (version == "1.1") {

    #if defaultTimeZone is not specified, then read it for each value
    if (N > bigData) { print("processing dateTimeUTC...") }

    if (is.null(time_diff)) {
      DateTimeUTC = xpathSApply(doc, "//sr:value", xmlGetAttr, name="dateTimeUTC", namespaces=ns)
      DateTimeUTC <- as.POSIXct(DateTimeUTC, format="%Y-%m-%dT%H:%M:%S", tz="GMT")
      UTCOffset = xpathSApply(doc, "//sr:value", xmlGetAttr, name="timeOffset", namespaces=ns)
      UTCOffset <- ifelse(grepl(":", UTCOffset),
                          as.numeric(substr(UTCOffset, nchar(UTCOffset)-4, nchar(UTCOffset)-3)),
                          as.numeric(UTCOffset))
      utcDiff = as.difftime(UTCOffset, units="hours")
      DateTime = as.POSIXct(DateTimeUTC + utcDiff)
      if (UTCOffset[1] > 0) {
        attr(DateTime, "tzone") <- paste("Etc/GMT+", UTCOffset[1], sep="")
      }
      if (UTCOffset[1] < 0) {
        attr(DateTime, "tzone") <- paste("Etc/GMT", UTCOffset[1], sep="")
      }
    } else {
      DateTime <- xpathSApply(doc, "//sr:value", xmlGetAttr, name="dateTime", namespaces=ns)
      zone="GMT"
      if (as.numeric(diff_text) > 0) {
        zone <- paste("Etc/GMT+", as.numeric(diff_text), sep="")
      }
      if (as.numeric(diff_text) < 0) {
        zone <- paste("Etc/GMT", as.numeric(diff_text), sep="")
      }
      DateTime <- as.POSIXct(DateTime, format="%Y-%m-%dT%H:%M:%S", tz=zone)
      UTCOffset = rep(as.numeric(diff_text), N)
      DateTimeUTC = DateTime - time_diff
      attr(DateTimeUTC, "tzone") <- "GMT"
    }

    if (N > bigData) { print("processing methodCode...") }
    methodCode = xpathSApply(doc, "//sr:value", xmlGetAttr, name="methodCode", namespaces=ns)
    methodCode <- unlist(methodCode)
    if (is.null(methodCode)) { methodCode <- NA }

    if (N > bigData) { print("processing sourceCode...") }
    sourceCode = xpathSApply(doc, "//sr:value", xmlGetAttr, name="sourceCode", namespaces=ns)
    sourceCode <- unlist(sourceCode)
    if (is.null(sourceCode)) { sourceCode <- NA }

    if (N > bigData) { print("processing qualityControlLevelCode...") }
    qcCode = xpathSApply(doc, "//sr:value", xmlGetAttr, name="qualityControlLevelCode", namespaces=ns)
    qcCode <- unlist(qcCode)
    if (is.null(qcCode)) { qcCode <- NA }

    nodata = as.numeric(xpathSApply(doc, "//sr:noDataValue", xmlValue, namespaces=ns))

  } else {

    #WaterML 1.0 usually doesn't provide information on UTC offset
    if (N > bigData) { print("processing dateTime...") }
    DateTimeUTC = xpathSApply(doc, "//sr:value", xmlGetAttr, name="dateTime", namespaces=ns)
    DateTime <- DateTimeUTC
    UTCOffset <- rep(0, N)

    if (N > bigData) { print ("processing methodID...")}
    methodCode <-  xpathSApply(doc, "//sr:value", xmlGetAttr, name="methodID", namespaces=ns)
    methodCode <- unlist(methodCode)
    if (is.null(methodCode)) { methodCode <- NA }

    if (N > bigData) { print ("processing sourceID...")}
    sourceCode <- xpathSApply(doc, "//sr:value", xmlGetAttr, name="sourceID", namespaces=ns)
    sourceCode <- unlist(sourceCode)
    if (is.null(sourceCode)) { sourceCode <- NA }

    if (N > bigData) { print ("processing qualityControlLevel...")}
    qcCode = xpathSApply(doc, "//sr:value", xmlGetAttr, name="qualityControlLevel", namespaces=ns)
    qcCode <- unlist(qcCode)
    if (is.null(qcCode)) { qcCode <- NA }
    if (length(qcCode) < N) { qcCode <- NA }

    nodata = as.numeric(xpathSApply(doc, "//sr:NoDataValue", xmlValue, namespaces=ns))
  }

  #make the data frame
  df <- data.frame(
    time=DateTime,
    DataValue=as.numeric(val),
    UTCOffset=UTCOffset,
    Qualifier=qualifier,
    CensorCode=censorCode,
    DateTimeUTC=DateTimeUTC,
    MethodCode=methodCode,
    SourceCode=sourceCode,
    QualityControlLevelCode=qcCode,
    stringsAsFactors=TRUE
  )

  if (nrow(df) == 0) {
    print("NOTE: No data values found in this time series")
    end.parse.time <- Sys.time()
    parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
    attr(df, "parse.status") <- "NO_VALUES_FOUND"
    attr(df, "parse.time") <- parse.time
    return(df)
  }

  #normal case: no aggregation
  df[df$DataValue == nodata,2] <- NA

  #special case: daily data aggregation
  if (!is.null(daily)) {
    validdata <- na.omit(df)
    if (nrow(validdata) == 0) {
      end.parse.time <- Sys.time()
      parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
      attr(df, "parse.status") <- "NO_VALUES_FOUND"
      attr(df, "parse.time") <- parse.time
      print("no valid data found!")
      return(df)
    }
    validdata$time <- as.Date(as.POSIXct(validdata$time))
    dailyValues <- aggregate(validdata$DataValue, list(validdata$time), daily)
    names(dailyValues)[1] <- "time"
    names(dailyValues)[2] <- "DataValue"

    end.parse.time <- Sys.time()
    parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
    attr(dailyValues, "download.status") <- attr(df, "download.status")
    attr(dailyValues, "download.time") <- attr(df, "download.time")
    attr(dailyValues, "parse.status") <- "OK"
    attr(dailyValues, "parse.time") <- parse.time
    return(dailyValues)
  }

  end.parse.time <- Sys.time()
  parse.time <- as.numeric(difftime(end.parse.time, begin.parse.time, units="sec"))
  attr(df, "download.time") <- download.time
  attr(df, "download.status") <- download.status
  attr(df, "parse.status") <- "OK"
  attr(df, "parse.time") <- parse.time
  return(df)
}
