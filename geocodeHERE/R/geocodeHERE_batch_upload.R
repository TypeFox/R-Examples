
#' Bulk upload addresses to Nokia HERE batch geocoding API
#'
#' Use the Nokia HERE batch geocoding API to geocode lots of addresses in a single call instead of looping over the addresses one by one.
#' @param address_string Character string containing the addresses to be geocoded. Output from 'format_vec_for_upload(...)
#' @param email_address Character string containing an email address. Nokia will email you here when the job is done.
#' @param App_id App_id to use the production HERE API. Get one here... http://developer.here.com/get-started. If left blank, will default to demo key with an unknown usage limit.
#' @param App_code App_code to use the production HERE API. Get one here... http://developer.here.com/get-started. If left blank, will default to demo key with an unknown usage limit.
#' @param quiet TRUE / FALSE indicating whether to write the POST information to the console
#' @return request_id as a string
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
geocodeHERE_batch_upload <- function(address_string, email_address, App_id="",
                                     App_code="", quiet=TRUE){
  if(!is.character(address_string)){stop("'address_string' must be a character string")}
  if(!is.character(email_address)){stop("'file_name' must be a character string")}
  if(!is.character(App_id)){stop("'App_id' must be a character string")}
  if(!is.character(App_code)){stop("'App_code' must be a character string")}

  if(App_id=="" & App_code==""){
    App_id <- "inwresuveWra5ebewaSweh"
    App_code <- "zBWCuMTr-PrXwr6pc5uqLg"    
    base_url <- "http://batch.geocoder.cit.api.here.com/6.2/jobs"
  }else{
    base_url <- "http://batch.geocoder.api.here.com/6.2/jobs"
  }

  v <- ifelse(quiet, httr::verbose(), NULL)
  a <- httr::POST(base_url, encode="multipart",
            body=address_string,
            query=list(
              action="run",
              mailto=email_address,
              maxresults="1",
              language="es-ES",
              header="true",
              indelim="|",
              outdelim="|",
              outcols="displayLatitude,displayLongitude,houseNumber,street,district,city,postalCode,county,state,country,matchLevel,relevance", # i shortened this for the example
              outputCombined="false",
              app_code=App_code,
              app_id=App_id),
            v)
  response <- httr::content(a)

  if(length(response$Response) > 0){
    request_id <- response$Response$MetaInfo$RequestId
    ret <- request_id
  }else{
    stop(paste("ERROR: ", response))
  }

  return(ret)
}
