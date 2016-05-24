ExtractData <-
function(sInputFile,
               sQuerySTARTTIME = NULL,sQueryENDTIME = NULL,
               sQueryTransmitterList = NULL,sQueryReceiverList = NULL,
               sQueryStationList = NULL, sQueryDataType = NULL)
{
  fEvaluateQuery <- sInputFile
  
  # query time
  if (is.null(sQuerySTARTTIME) == FALSE)
  {
    s_dates <- as.POSIXct(strptime(c(sQuerySTARTTIME,sQueryENDTIME), format = "%Y-%m-%d %H:%M:%S"))
    fEvaluateQuery <- fEvaluateQuery[fEvaluateQuery[,1] %in% s_dates[1]:s_dates[2],]
  }

  # query transmitter
  if (is.null(sQueryTransmitterList) == FALSE)
    fEvaluateQuery <- fEvaluateQuery[fEvaluateQuery[,2] %in% sQueryTransmitterList,]
 
  # query receiver
  if (is.null(sQueryReceiverList) == FALSE)
    fEvaluateQuery <- fEvaluateQuery[fEvaluateQuery[,5] %in% sQueryReceiverList,]

  # query station
  if (is.null(sQueryStationList) == FALSE)
    fEvaluateQuery <- fEvaluateQuery[fEvaluateQuery[,6] %in% sQueryStationList,]

  # query data type
  if (is.null(sQueryDataType) == FALSE)
    fEvaluateQuery <- fEvaluateQuery[fEvaluateQuery[,4] %in% sQueryDataType,]
  
  return(fEvaluateQuery)
}
