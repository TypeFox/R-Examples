mixpanelGetEventsFromFiles = function(
  account,
  dateFrom,
  dateTo=dateFrom,
  eventNames=c(),       # Import only selected events?
  select=TRUE,          # Requested column names.
  blocksize=500000,     # Block size in reading files from disk. 
  verbose=TRUE
) {
  dates = as.character(seq(as.Date(dateFrom), as.Date(dateTo), "day"))
  alldata = matrix(NA, 0, 0)
  
  for (date in dates) {
    cat("#########", date, "########", date(), "\n")
    
    txtFilePath = paste(account$dataPath, "/", "Events-", date, ".txt", sep="")
    con <- file(txtFilePath, "r", blocking=TRUE)
    
    block = 1
    while (TRUE) {
      
      cat("### Scan block", block, "of events  -", date(), "\n")
      data = readLines(con, n=blocksize)
      if(length(data) == 0)
        break
      
      cat("### Subset, parse & merge   -", date(), "\n")
      if (length(eventNames) > 0) {
        inds = rep(FALSE, length(data))
        for (name in eventNames) 
          inds =  inds | grepl(name, data, fixed=TRUE)
        data = data[inds]
      }
      
      newdata = eventsJson2RMatrix(data, select)
      
      if (length(eventNames) > 0) {
        ## Only needed in the case event names were not handled by grep.
        ind = (newdata[, "event"] %in% eventNames)
        alldata = merge.matrix(alldata, newdata[ind, , drop=FALSE])
      } else {
        alldata = merge.matrix(alldata, newdata)
      }
      
      block = block + 1
    }
    close(con)
  }
  
  cat("### Flatten matrix           -", date(), "\n")
  alldata = getFlatMatrix(alldata)
  cat("### Done.                    -", date(), "\n")
  alldata
}
