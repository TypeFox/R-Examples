
#' Gets the results of a batch geocoding job
#'
#' Download the result to a temp file, extract zip, read pipe delimited file, return as dataframe
#' @param request_id Character string containing a request_id. This is returned from geocodeHERE_batch_upload(...)
#' @param App_id App_id to use the production HERE API. Get one here... http://developer.here.com/get-started. If left blank, will default to demo key with an unknown usage limit.
#' @param App_code App_code to use the production HERE API. Get one here... http://developer.here.com/get-started. If left blank, will default to demo key with an unknown usage limit.
#' @return Dataframe of results from the geocoding job
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
geocodeHERE_batch_get_data <- function(request_id="", App_id="", App_code=""){
  if(!is.character(request_id)){stop("'request_id' must be a character string")}
  if(request_id==""){stop("'request_id' must be have a value")}
  if(geocodeHERE_batch_status(request_id) != "completed"){
    stop("Batch geocoding is not completed yet")}

  if(App_id=="" & App_code==""){
    App_id <- "inwresuveWra5ebewaSweh"
    App_code <- "zBWCuMTr-PrXwr6pc5uqLg"    
    base_url <- "http://batch.geocoder.cit.api.here.com/6.2/jobs"
  }else{
    base_url <- "http://batch.geocoder.api.here.com/6.2/jobs"
  }

  download_url <- paste0(base_url, "/",
                         request_id,
                         "/result",
                         "?app_id=", App_id,
                         "&app_code=", App_code)
  file_path <- tempfile()
  a <- httr::GET(download_url, httr::write_disk(file_path, overwrite=TRUE))
  response <- httr::content(a)

  if(is.list(response)){
    stop(paste("ERROR: ", response$Details))
  }

  con <- utils::unzip(file_path, exdir=paste0(file_path, "tmp"))
  df <- utils::read.delim(con, stringsAsFactors=F, sep="|")
  unlink(file_path)

  return(df)
}
