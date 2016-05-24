gaRequest <-
function(id, dimensions=NA, metrics, start, end, token=NA, sort=NA, max=1000, segment=NA, filters=NA, allResults=FALSE) {
   
  ##INTERNAL FUNCTION BUILD AND FETCH --------------------------------------------------------------------------------------------------------
  buildAndFetch <- function(q){
    ##BUILD QUERY AND FETCH RESULTS
    r <-list()
    for (x in 1:length(q)) {
      r[[x]] <- GET("https://www.googleapis.com/analytics/v3/data/ga",
                    query = q[[x]]
      )
    }
    
    #PARSE JSON RESPONSE FROM GA
    for (y in 1:length(r)) {
      r[[y]] <- jsonlite::toJSON(content(r[[y]]))
      r[[y]] <- jsonlite::fromJSON(r[[y]], flatten=TRUE)
    }
    return(r)
  }
  
  ##INTERNAL FUNCTION LIST TO DF --------------------------------------------------------------------------------------------------------
  toDF <- function(r) {
    
    df <- list()
    for (x in 1:length(r)) {
      df[x] <- list(cbind(
        data.frame(r[[x]]$profileInfo, stringsAsFactors=FALSE),
        data.frame(r[[x]]$query['start-date'], stringsAsFactors=FALSE),
        data.frame(r[[x]]$query['end-date'], stringsAsFactors=FALSE),
        data.frame(r[[x]]$totalResults, stringsAsFactors=FALSE),
        data.frame(r[[x]]$rows, stringsAsFactors=FALSE)
      ))
    }
    
    df <- do.call('rbind', df)
    
    ##RENAME COLUMNS
    colnames(df) <- c(
      names(r[[1]]$profileInfo),
      names(r[[1]]$query['start-date']),
      names(r[[1]]$query['end-date']),
      'totalResults',
      unlist(r[[1]]$columnHeaders$name)
    )
    
    #FORMAT METRIC DATA AS NUMERIC
    metricNames <- as.vector(r[[1]]$query$metrics)
    for (i in metricNames){
      df[,i] <- as.numeric(df[,i])
    }
    
    ##REMOVE GA FROM COLUMN NAMES
    colnames(df) <- gsub('ga:','',colnames(df))
    
    return(df)
  }
  
  
  ##PRIORITIZES SUPPLIED TOKEN FIRST. THEN CHECKS FOR ENV VARIABLE, THEN DEFAULTS TO NA. ------------------------------------------------------------
  if (!is.na(token)){
    token <- token
  } else if (exists('GAR_ACCESS_TOKEN', envir=envGAR)){
    token <- get('GAR_ACCESS_TOKEN', envir=envGAR)
  } else {
    token <- token
  }
  
  
  ##CREATE LIST OF INITIAL QUERY PARAMETERS --------------------------------------------------------------------------------------------------------
  queryList <- as.list(id)
  for (x in 1:length(id)) { 
    queryList[[x]] <- list(
      'ids'=id[x],
      'dimensions'=dimensions,
      'metrics'=metrics,
      'start-date'=start,
      'end-date'=end,
      'sort'=sort,
      'max-results'=max,
      'filters'=filters,
      'segment'=segment,
      'access_token'=token
    )
  }
    
  ##REMOVE ANY UNUSED PARAMETERS
  for (x in 1:length(id)){
    queryList[[x]] <- queryList[[x]][!is.na(queryList[[x]])]
  }
    
  
  ##BUILD AND FETCH INITIAL QUERY --------------------------------------------------------------------------------------------------------
  finalDf <- buildAndFetch(queryList)
  
  
  ###CHECK FOR INITIAL ERRORS
  if(TRUE %in% grepl('error', names(finalDf[[1]]) )){
    ##IF ERROR THEN DATA FRAME OF ERROR RESPONSE CODE AND MESSAGE
    
    tmp <- list()
    for (x in 1: length(id)){
      tmp[x] <- list(data.frame(finalDf[[x]]))
    }
    
    finalDf <- do.call('rbind',tmp)
    
  } else {
  
  finalDf <- toDF(finalDf)
  } #END IF ELSE STATEMENT
  
  
  ##CREATE LIST OF PAGINATION QUERY PARAMETERS --------------------------------------------------------------------------------------------------------
  if (allResults==TRUE & max==10000) {
  
  totalResults <- aggregate(data=finalDf, totalResults~tableId, FUN=mean)
  totalResults$pages <- floor(totalResults$totalResults/max)
  
  tmp <- data.frame(profileId=NA, pages=NA)
  for (x in totalResults$tableId){
    tmp <- rbind(tmp,
            data.frame(profileId=totalResults$tableId[totalResults$tableId==x], pages=0:totalResults$pages[totalResults$tableId==x], stringsAsFactors=FALSE, row.names=NULL)
    )
  }

  tmp <- tmp[!is.na(tmp$pages) & tmp$pages>0,]
  rm(totalResults)
  
  queryListPaginate <- list()
  for (x in 1:length(tmp$profileId)) {
    queryListPaginate[[x]] <- list(
      'ids'=tmp$profileId[x],
      'dimensions'=dimensions,
      'metrics'=metrics,
      'start-date'=start,
      'end-date'=end,
      'sort'=sort,
      'max-results'=max,
      'filters'=filters,
      'segment'=segment,
      'start-index'=(tmp$pages[x]*max)+1,
      'access_token'=token
    )
  }
  
  ##REMOVE ANY UNUSED PARAMETERS
  for (x in 1:length(queryListPaginate)){
    queryListPaginate[[x]] <- queryListPaginate[[x]][!is.na(queryListPaginate[[x]])]
  }
  
  ##BUILD AND FETCH INITIAL QUERY --------------------------------------------------------------------------------------------------------
  paginateDf <- buildAndFetch(queryListPaginate)
  paginateDf <- toDF(paginateDf)
  
  finalDf <- rbind(finalDf, paginateDf)
  } else {
  finalDf <- finalDf
  }
  
  ##REMOVE UNNECESSARY COLUMNS FOR FINAL DF ------------------------------------------------------------------------------------------------
  finalDf <- finalDf[,colnames(finalDf) != 'totalResults']
  return(finalDf)
  }