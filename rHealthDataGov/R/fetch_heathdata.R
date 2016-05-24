# Retrieve data from the healthdata.gov data API.
# More info here: http://www.healthdata.gov/data-api


.construct_json_string <- function(resource_id, filter, offset=0){
  # This constructs the initial json string for a query with a given filter set
  
  if (!is.null(filter)){
    filter_vec <- mapply(names(filter), filter, FUN=function(k,v) sprintf('"%s": "%s"', k, v))
    filter_str <- sprintf('{%s}', paste(filter_vec, collapse=","))
    jsontext <- sprintf('{"resource_id": "%s","filters": %s, "offset": %s}', resource_id, filter_str, offset)
  } else {
    jsontext <- sprintf('{"resource_id": "%s", "offset": %s}', resource_id, offset)  #no filter
  }
  return(jsontext)
}

.quick_fetch <- function(jsontext){
  # This will fetch results, but will be limited by the API at max = 100 results
  
  api_url <- "http://hub.Healthdata.gov/api/action/datastore_search"
  req <- POST(api_url, body = jsontext, add_headers("Content-type" = "application/json"))
  
  if (req$status_code >= 400) {
    stop(sprintf("HealthData.gov API returned an error: HTTP status code %s, %s", req$status_code, req$headers$statusmessage))
  }
  stop_for_status(req) 
  reslist <- content(req, "parsed")
  
  if (!reslist$success) stop("HealthData.gov API returned an error.")
  if (length(reslist$result$records)==0) stop("Zero records match your filter. Nothing to return.\n See the 'filters' data object from this package for valid filter values.")
  
  return(reslist)
}


.fetch_results <- function(jsontext, resource_id, filter){
  # Fetch the first batch of results (limit 100 records)
  
  reslist <- .quick_fetch(jsontext)
  total_records <- reslist$result$total  #Total number of records that match the query
  records <- reslist$result$records  #Records returned by API call
  field_names <- sapply(reslist$result$fields, function(x) x$id)  #Field names
  field_types <- sapply(reslist$result$fields, function(x) x$type)  #Corresponding field types
  
  # If you have not yet retrieved all records, calculate the # of remaining calls required
  extra_calls <- ifelse((length(records) < total_records), floor(total_records/length(records)), 0)
  
  if (extra_calls>0){
    # Might update the following to an apply or parSapply function (with a new parallel=TRUE arg)
    all_records <- list(records)
    for (i in seq(extra_calls)) {
      # Keep making API requests with an increasing offset value until you get all the records
      # Append new records to existing `records` list
      api_hardlimit <- 100  #healthdata.gov has a hard limit of 100 records per request
      jsontext <- .construct_json_string(resource_id=resource_id, filter=filter, offset=api_hardlimit*i)
      all_records[[i+1]] <- .quick_fetch(jsontext)$result$records
    }
    records <- unlist(all_records, recursive=FALSE)
  }
  
  # Replace NULL with NA so we can convert into a data.frame    
  records <- lapply(records, function(x) lapply(x, function(x) ifelse(is.null(x), NA, x)))     
  df <- data.frame(matrix(unlist(records), nrow=length(records), byrow=T), stringsAsFactors=FALSE)
  names(df) <- names(records[[1]])
  
  # Convert column types
  for (col in names(df)){
    field_type <- field_types[which(field_names==col)]
    if (field_type=="numeric"){
      df[[col]] <- as.numeric(df[[col]])
    } else if (field_type=="int4"){
      df[[col]] <- as.integer(df[[col]])
    } else if (field_type=="int8"){
      df[[col]] <- as.integer64(df[[col]])
    } else if (field_type=="timestamp"){
      df[[col]] <- as.POSIXct(df[[col]], tz="UTC", format="%Y-%m-%dT%H:%M:%S")
    } else if (field_type!="text"){
      print(sprintf("Unknown field type: %s", field_type))
    }
  }  
  
  return(df)
}




fetch_healthdata <- function(resource="hosp", filter=NULL){
    # Return a data.frame of the records that match your query
    
    resource_id <- resources$resource_id[which(resources$resource==resource)]
    if (length(resource_id)==0) {
      stop("This resource does not exist.  Please input a character string that is one of the resources$resource labels.")
    }
    if (!is.null(filter)) {
      bad_filter_names <- sapply(names(filter), function(x) !(x %in% names(filters[[resource]])))
      if (sum(bad_filter_names)>0) {
        baddies <- paste(names(filter)[bad_filter_names], collapse=", ")
        stop(paste("The following filter name(s) do not exist for this resource", baddies, sep=": "))  
      }      
    }
    jsontext <- .construct_json_string(resource_id=resource_id, filter=filter, offset=0)
    df <- .fetch_results(jsontext=jsontext, resource_id=resource_id, filter=filter) 
    return(df)
}
