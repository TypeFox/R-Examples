mixpanelGetProfiles <- function(
  account,
  where='',
  select=TRUE,       # Select only these columns. <TRUE> selects all columns.
  maxPage=100000,    # Import until maxPage. Stop when a page does not exist. 
  
  verbose=TRUE
) {
  args = list(where=where)
  
  if(verbose)
    cat("*** Page 0 ***\n")
  data = mixpanelGetData(account, "engage/", args, data=TRUE, verbose=verbose)
  data = jsonlite::fromJSON(data)
  
  args$session_id = data$session_id
  pageSize = data$page_size
  totalCount = data$total
  if (totalCount == 0) {
    stop("No profile found.")
  }
  if(verbose)
    cat("*** Total profiles count is", totalCount, "\n")
  
  alldata = profilesJson2RMatrix(data, select)
  
  for (page in 1:maxPage) {
    if(verbose)
      cat("*** Page", page, "***\n")
    args$page = page
    if (page * pageSize >= totalCount || maxPage <= 0) break
    data = mixpanelGetData(account, "engage/", args, data=TRUE, verbose=verbose)
    data = jsonlite::fromJSON(data)
    newdata = profilesJson2RMatrix(data, select)
    alldata = merge.matrix(alldata, newdata)
  }
  if(verbose)
    cat("*** Finished.\n")
  
  alldata
}
