
#' @title Data Item Explorer API
#' @description This function connects to the UK National Grid's API for Data Item Explorer, which is a major data source for gas-related information. Internet connection must be available.
#' @details The function submits an enquiry to the destination API using XML over HTTP using SOAP standard. The HTTP response is in XML format which function will parse internally and return a R dataframe object.
#' @param dataitems A vector of characters containing data items to enquire via API.
#' @param fromdate A character object specifying the start date. Date is inclusive.
#' @param todate A character object specifying the end date. Date is inclusive.
#' @param datetype A character object specifying the data type. Defaults to \code{gasday}
#' @param latestflag A character object with length of one to specify whether to extract the latest data. This can either be \code{Y} or \code{N}. Defaults to \code{Y}.
#' @param applicableforflag A character object with length of one to specify whether dates specified are 'applicable for' or 'applicable on'. This can either be \code{Y} or \code{N} where \code{Y} indicates 'applicable for'. Defaults to \code{Y}.
#' @param apiurl A character object which points to National Grid's SOAP API. Under most circumstances users do not have to change this. Defaults to 'http://marketinformation.natgrid.co.uk/MIPIws-public/public/publicwebservice.asmx'
#' @return A dataframe object containing API response data.
#' @examples
#' # Specify the data item(s) to enquire from API
#' dataitems <- c('Storage Injection, Actual',
#'                'Storage Withdrawal, Actual')
#'
#' # Initialise API (requires internet connection for this step)
#' response <- dataItemExplorer(dataitems,
#'                              fromdate = '2013-10-01',
#'                              todate='2015-09-30')
#'
#' # Visualise the results on a chart
#' library(ggplot2)
#' ggplot(response,aes(x=ApplicableFor,y=Value,colour=PublicationObjectName)) + geom_line()
#' @author Timothy Wong, \email{timothy.wong@@hotmail.co.uk}
#' @references
#' \itemize{
#' \item Graphical User Interface for Data Item Explorer\cr
#' \url{http://www2.nationalgrid.com/uk/industry-information/gas-transmission-operational-data/data-item-explorer/}
#' \item API specification\cr
#' \url{http://marketinformation.natgrid.co.uk/MIPIws-public/public/publicwebservice.asmx?op=GetPublicationDataWM}
#' }
#' @export
dataItemExplorer<- function(dataitems,
                            fromdate,
                            todate,
                            datetype='gasday',
                            latestflag='Y',
                            applicableforflag='Y',
                            apiurl = 'http://marketinformation.natgrid.co.uk/MIPIws-public/public/publicwebservice.asmx') {

  # Creates SOAP XML request
  soap.request <- base::paste0('<?xml version="1.0" encoding="utf-8"?><soap12:Envelope xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:soap12="http://www.w3.org/2003/05/soap-envelope"><soap12:Body><GetPublicationDataWM xmlns="http://www.NationalGrid.com/MIPI/"><reqObject><LatestFlag>',latestflag,'</LatestFlag><ApplicableForFlag>',applicableforflag,'</ApplicableForFlag><ToDate>',todate,'</ToDate><FromDate>',fromdate,'</FromDate><DateType>',datetype,'</DateType><PublicationObjectNameList>',paste0('<string>',dataitems,'</string>', collapse = ''),'</PublicationObjectNameList></reqObject></GetPublicationDataWM></soap12:Body></soap12:Envelope>')

  # Initialises a text gatherer
  h <- RCurl::basicTextGatherer()
  h$reset()

  # Initiates SOAP request via HTTP
  RCurl::curlPerform(url = apiurl, httpheader = c(Accept="text/xml", 'Content-Type'='application/soap+xml; charset=utf-8'),postfields = soap.request, verbose=TRUE, writefunction = h$update)

  # Writes SOAP response into character
  soap.response <- h$value()

  # converts SOAP response into an XML object
  soap.response.doc <- XML::xmlTreeParse(soap.response, replaceEntities = TRUE , useInternalNodes = TRUE)

  # converts SOAP response into a list of objects
  objects <- (XML::xmlToList(soap.response.doc,simplify = TRUE))[[1]][[1]]

  list.of.data.frames<-list()

  for(i in 1:length(objects)) {
    # Porcessing data frames
    rows <- objects[[i]][['PublicationObjectData']]

    list.of.rows <- base::data.frame()

    for(r in 1:length(rows)){
      # Processing individual row
      list.of.rows <- base::rbind(list.of.rows,as.data.frame(rows[r][[1]],stringsAsFactors = FALSE))
    }
    list.of.rows$PublicationObjectName <- objects[[i]][[1]]

    list.of.data.frames <- base::rbind(list.of.data.frames, as.data.frame(list.of.rows))
  }

  # Convert all columns to appropiate types before returning data frame
  result <- base::within(data = list.of.data.frames,expr = {
    ApplicableAt <- strptime(ApplicableAt,'%Y-%m-%dT%H:%M:%SZ','UTC')
    ApplicableFor <- strptime(ApplicableFor,'%Y-%m-%dT%H:%M:%SZ','UTC')
    Value <- as.numeric(Value)
    GeneratedTimeStamp <- strptime(GeneratedTimeStamp,'%Y-%m-%dT%H:%M:%SZ','UTC')
    QualityIndicator <- factor(QualityIndicator)
    Substituted <- factor(QualityIndicator)
    CreatedDate <- strptime(CreatedDate,'%Y-%m-%dT%H:%M:%SZ','UTC')
    PublicationObjectName <- as.factor(PublicationObjectName)
  })

  return (result)
}
