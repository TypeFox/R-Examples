
# Open new chart's window using this chartId
#
getRemoteChartUrl <-
  function(chartId='',openBrowser=TRUE)
  {
    chartUrl<- buildURL(paste( "/remoteChart?username=",machina$session$userId,sep = ""))
    if(chartId != '')
    {
      chartUrl<- paste(chartUrl,"?chartId=",chartId,sep="")
    }
    if(openBrowser)
    {
      browseURL(chartUrl)
    }
    
    return(chartUrl)
  }
# Opens the browser only if it's not already opened
openRemoteChartForFirstTime<-
  function(chartId='',verbose = FALSE)
  {
    #remoteChartConnected
    
    rsp <-
      POST(
        buildURL("/api/user/remoteChart/remoteChartConnected"),
        add_headers(access_token = machina$session$accessToken),
        body = list(chartId=chartId),
        encode="json",
        if (verbose) { verbose() })
    
    if (rsp$status_code != 200)
    {
      stop("[", rsp$status_code, "]: ", rawToChar(rsp$content))
    }
    result=jsonlite::fromJSON(rawToChar(rsp$content))
    if(!result$result)
    {
     return(getRemoteChartUrl(openBrowser = TRUE))
    }
  }

# Draw the days to a remote chart
#
#   Parameters:
#     chartId - The unique ID that connect the webpage using socke.io
#     seriesDataList - The seriesDataList, for example list(GOOG=c(...),SPY=c(...))
drawSeries <-
  function(seriesDataList,chartId='', verbose = FALSE)
  {
    rsp <-
      POST(
        buildURL("/api/user/remoteChart/sendMessagToRemoteChart"),
        add_headers(access_token = machina$session$accessToken),
        body = list(chartId=chartId,
                    seriesData=seriesDataList,
                    action='drawSeries'),
        encode="json",
        if (verbose) { verbose() })
    
    if (rsp$status_code != 200)
    {
      stop("[", rsp$status_code, "]: ", rawToChar(rsp$content))
    }
    openRemoteChartForFirstTime();
  }

# Draw these rows to a remote chart
#
#   Parameters:
#     chartId - The unique ID that connect the webpage using socke.io
#     rowsArray - The rows indexes
drawRows <-
  function(rowsArray,startDate = NULL,endDate = NULL,chartId='', verbose = FALSE)
  {
    seriesList<- list();
    for(item in rowsArray){
      rowData<- getRow(item,startDate = startDate,endDate = endDate,includeData = TRUE,verbose = verbose)
      #chartKey<- paste('Chart',item)
      chartKey<- machina$session$model$query[item]
      seriesList[[chartKey]]<- rowData$days
    }
    drawSeries(seriesList)
  }