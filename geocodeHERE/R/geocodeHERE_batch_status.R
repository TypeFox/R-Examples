
#' Check the status on a batch geocoding job
#'
#' Check the status on a batch geocoding job
#' @param request_id Character string containing a request_id. This is returned from geocodeHERE_batch_upload(...)
#' @param full_list TRUE / FALSE indicating whether to return the full response from Nokia HERE or just the "Status" portion of the response
#' @param App_id App_id to use the production HERE API. Get one here... http://developer.here.com/get-started. If left blank, will default to demo key with an unknown usage limit.
#' @param App_code App_code to use the production HERE API. Get one here... http://developer.here.com/get-started. If left blank, will default to demo key with an unknown usage limit.
#' @return A string representing the status of the job. The full reply from Nokia as a list if the full_list parameter is TRUE
#' @keywords geocode batch
#' @export
#' @examples
#' addresses <- chicago_landmarks[,"Address"]
#' addresses <- paste(addresses, "chicago IL")
#' addresses_df <- data.frame(id=1:length(addresses), addresses=addresses)
#' address_str <- df_to_string(addresses_df)
#' \dontrun{
#' request_id <- geocodeHERE_batch_upload(address_string = address_str,
#'                                        email_address = "youremail<at>domain.com")
#' geocodeHERE_batch_status(request_id)
#' geocode_data <- geocodeHERE_batch_get_data(request_id)
#' addresses_df <- merge(addresses_df, geocode_data, by.x="id", by.y="recId", all.x=T)
#' }
geocodeHERE_batch_status <- function(request_id="", full_list=FALSE, App_id="", App_code=""){
  if(!is.character(request_id)){stop("'request_id' must be a character string")}
  if(request_id==""){stop("'request_id' must be have a value")}
  if(!is.logical(full_list)){stop("'full_list' must be a logical value")}

  if(App_id=="" & App_code==""){
    App_id <- "inwresuveWra5ebewaSweh"
    App_code <- "zBWCuMTr-PrXwr6pc5uqLg"    
    base_url <- "http://batch.geocoder.cit.api.here.com/6.2/jobs"
  }else{
    base_url <- "http://batch.geocoder.api.here.com/6.2/jobs"
  }

  status_url <- paste0(base_url, "/",
                       request_id,
                       "?action=status",
                       "&app_id=", App_id,
                       "&app_code=", App_code)
  a <- httr::GET(status_url)
  response <- httr::content(a)

  if(length(response$Response) > 0){
    status <- response$Response$Status
    if(full_list){
      ret <- response
    }else{
      ret <- status
    }
  }else{
    stop(paste("ERROR: ", response))
  }

  return(ret)
}
