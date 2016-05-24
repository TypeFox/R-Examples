#' HISCentral_GetSeriesCatalog
#'
#' This function searches the table of time series from the HIS Central catalog
#'
#' @import XML
#' @import httr
#' @param west The west longitude of the geographic
#'  bounding box in decimal degrees. Allowed values are between -180.0 and +180.0
#' @param south The south latitude of the geographic
#'  bounding box in decimal degrees. Allowed values are between -90.0 and +90.0
#' @param east The east longitude of the geographic
#'  bounding box in decimal degrees. Allowed values are between -180.0 and +180.0
#' @param north The north latitude of the geographic
#'  bounding box in decimal degrees. Allowed values are between -90.0 and +90.0
#' @param serviceID (optional): The ID of the service on HIS Central. To get the service ID,
#'  use the id column in the output of the GetServices() function.
#' @param keyword (optional): The concept keyword (common name of variable) for
#'  searching the sites on HIS Central. Examples include Temperature, Precipitation, Snow Depth,... If the Keyword is not
#'  specified then sites with any variable will be returned.
#' @param beginDate (optional): The begin date of the observations in yyyy-mm-dd format.
#' @param endDate (optional): The end date of the observations in yyyy-mm-dd format.
#' @return a data.frame of series catalog entries. The data.frame has the following columns:
#' \itemize{
#' \item ServiceCode: The code of the HydroServer
#' \item ServiceURL: The URL of the server. Use this as the server parameter in GetValues() function.
#' \item FullSiteCode: The complete unique code of the site in the format NETWORK:CODE.
#'               Use this value as the siteCode parameter in the GetValues function.
#' \item FullVariableCode: The complete unique code of the site in the format VOCABULARY:CODE.
#'               Use this value as the variableCode parameter in the GetValues function.
#' \item BeginDateTime: The local date/time of the first observation of the time series in POSIXct format.
#' \item EndDateTime: The local date/time of the last observation of the time series in POSIXct format.
#' \item ValueCount: The number of measurements in the time series
#' \item SiteName: The name of the site.
#' \item Latitude:  The WGS84 latitude of the site in decimal degrees
#' \item Longitude: The WGS84 longitude of the site in decimal degrees
#' \item DataType: The data type of the variable
#' \item ValueType: The type of the observation (field observation, sample, or derived value)
#' \item SampleMedium: The sample medium (air, water or other)
#' \item TimeUnits: The time units
#' \item TimeSupport: The length of the time period of one measurement
#' }
#' @keywords waterml
#' @export
#' @examples
#' #Getting all time series from the (14.1E, 49.9N, 14.3E, 50.1N) bounding box
#' series_catalog <- HISCentral_GetSeriesCatalog(west=14.1, south=49.9, east=14.3, north=50.1)

HISCentral_GetSeriesCatalog <- function(west, south, east, north,
                                serviceID=NULL, keyword=NULL, beginDate=NULL, endDate=NULL) {

  catalog = "http://hiscentral.cuahsi.org/webservices/hiscentral.asmx/GetSeriesCatalogForBox2"

  #create the URL
  servID = serviceID
  if (is.null(serviceID)) {
    servID=""
  }
  if (is.null(keyword)) {
    keyword=""
  }
  if (is.null(beginDate)) {
    beginDate="1900-01-01"
  }
  if (is.null(endDate)) {
    endDate="2050-01-01"
  }
  queryParameters <- list(xmin=west, ymin=south, xmax=east, ymax=north,
                          networkIDs=servID, conceptKeyword=keyword,
                          beginDate=beginDate, endDate=endDate)
  url <- paste(catalog, "?", "&xmin=", west, "&ymin=", south, "&xmax=", east,
               "&ymax=", north, "&networkIDs=", servID, "&conceptKeyword=", keyword,
               "&beginDate=", beginDate, "&endDate=", endDate,
               sep="")

  print(paste("searching sites from:", url, "..."))

  download.time <- system.time(
    tryCatch({
      downloaded <- FALSE
      response <- GET(catalog, query=queryParameters)
      downloaded <- TRUE
    },error=function(e){
      print(conditionMessage(e))
    })
  )

  if (!downloaded) {
    return(NULL)
  }

  status.code <- http_status(response)$category
  print(paste("download time:", download.time["elapsed"], "seconds, status:", status.code))

  ######################################################
  # Parsing the WaterML XML Data                       #
  ######################################################

  print("reading sites XML data...")
  doc <- tryCatch({
    xmlParse(response)
  }, warning = function(w) {
    print("Error reading HIS Central Data: Bad XML format.")
    return(NULL)
  }, error = function(e) {
    print("Error reading HIS Central Data: Bad XML format.")
    return(NULL)
  }
  )
  if (is.null(doc)) {
    return(NULL)
  }

  # specify the namespace information for HIS Central
  ns <- c(xsd="http://www.w3.org/2001/XMLSchema",
          xsi="http://www.w3.org/2001/XMLSchema-instance",
          sr="http://hiscentral.cuahsi.org/20100205/")

  # extract the data columns with XPath
  ServiceCode <- xpathSApply(doc, "//sr:ServCode", xmlValue, namespaces=ns)
  ServiceURL <- xpathSApply(doc, "//sr:ServURL", xmlValue, namespaces=ns)
  FullSiteCode <- xpathSApply(doc, "//sr:location", xmlValue, namespaces=ns)
  FullVariableCode <- xpathSApply(doc, "//sr:VarCode", xmlValue, namespaces=ns)
  BeginDateTime <- xpathSApply(doc, "//sr:beginDate", xmlValue, namespaces=ns)
  EndDateTime <- xpathSApply(doc, "//sr:endDate", xmlValue, namespaces=ns)
  SiteName <- xpathSApply(doc, "//sr:Sitename", xmlValue, namespaces=ns)
  Latitude <- xpathSApply(doc, "//sr:latitude", xmlValue, namespaces=ns)
  Longitude <- xpathSApply(doc, "//sr:longitude", xmlValue, namespaces=ns)
  DataType <- xpathSApply(doc, "//sr:datatype", xmlValue, namespaces=ns)
  ValueType <- xpathSApply(doc, "//sr:valuetype", xmlValue, namespaces=ns)
  SampleMedium <- xpathSApply(doc, "//sr:samplemedium", xmlValue, namespaces=ns)
  TimeUnits <- xpathSApply(doc, "//sr:timeunits", xmlValue, namespaces=ns)
  TimeSupport <- xpathSApply(doc, "//sr:TimeSupport", xmlValue, namespaces=ns)

  df <- data.frame(
    ServiceCode=ServiceCode,
    ServiceURL=ServiceURL,
    FullSiteCode=FullSiteCode,
    FullVariableCode=FullVariableCode,
    BeginDateTime=as.POSIXct(strptime(BeginDateTime, format="%m/%d/%Y")),
    EndDateTime=as.POSIXct(strptime(EndDateTime, format="%m/%d/%Y")),
    SiteName=SiteName,
    Latitude=as.numeric(Latitude),
    Longitude=as.numeric(Longitude),
    DataType=DataType,
    ValueType=ValueType,
    SampleMedium=SampleMedium,
    TimeUnits=TimeUnits,
    TimeSupport=as.numeric(TimeSupport),
    stringsAsFactors=FALSE
  )

  return(df)
}
