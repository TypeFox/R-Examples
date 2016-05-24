mixpanelGetEventsForProfiles <- function(
  account,
  distinctIDs=c("F05AD5E5-C4...", "437dbad6-72..."),    # Array of IDs.
  from="2010-01-01",   # ! Month numbers start HERE with 0!
  to="2020-01-01",     # ! Month numbers start HERE with 0!
  verbose=TRUE,
  ...     # Additional parameters, e.g. limit=5, ...
) {
## MP, 2015
##
  args = list(...)
  args$distinct_ids = arrayRtoJSON(distinctIDs)
  args$from_date = from
  args$to_date = to
  
  res = mixpanelGetData(account, "stream/query", args, data=TRUE, verbose=verbose)
  res = jsonlite::fromJSON(res)
  
  if ("status" %in% names(res) && res$status == "ok") {
    data = res$results$events
    if (length(data) > 0) {
      props = data[, 2]
      indID = which(colnames(props) == "distinct_id")
      data = cbind(distinctID=props[, indID], 
                   event=data[, 1], 
                   props[, -indID, drop=FALSE], 
                   stringsAsFactors=FALSE)
      data
    
    } else {
      cat("Event list empty")
      data.frame()
    }
    
  } else {
    print(res)
    data.frame()
  }
}
