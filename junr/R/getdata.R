#' @import httr
#' @import jsonlite
#'
NULL

#' Get a list of GUID's
#'
#' Get a list of all the available GUID's with datasets or views from the
#' base URL
#'
#' @param base_url The base URL of the Junar service
#' @param api_key The user's API key for the Junar service
#' @keywords GUID
#' @export
get_index <- function(base_url, api_key){
  if (missing(base_url)) {
    warning("Please add a valid base URL")
  }
  if (missing(api_key)) {
    warning("Please add a valid API key for the base URL you are trying to access")
  }
  try({
    r <- GET(paste(base_url, "?auth_key=", api_key, sep=""), accept_json())
    content_index <- fromJSON(content(r, "text"))
    return(content_index)
  })
}

#' Show GUID list
#'
#' Show a list of all available GUID's identifying the available data sets. It
#' only shows GUID's so that the list will fit in the console window.
#'
#' @param base_url The base URL of the Junar service
#' @param api_key The user's API key for the Junar service
#' @keywords GUID
#' @export
list_guid <- function(base_url, api_key){
  if (missing(base_url)) {
    warning("Please add a valid base URL")
  }
  if (missing(api_key)) {
    warning("Please add a valid API key for the base URL you are trying to access")
  }
  try({
    content_index <- get_index(base_url, api_key)
    return(content_index$guid)
  })
}

#' Show GUID titles
#'
#' Show a list of all available GUID titles. It only shows the titles so that
#' the list will fit in the console window.
#'
#' @param base_url The base URL of the Junar service
#' @param api_key The user's API key for the Junar service
#' @keywords GUID
#' @export
list_titles <- function(base_url, api_key){
  if (missing(base_url)) {
    warning("Please add a valid base URL")
  }
  if (missing(api_key)) {
    warning("Please add a valid API key for the base URL you are trying to access")
  }
  try({
    content_index <- get_index(base_url, api_key)
    return(content_index$title)
  })
}

#' Get data for a given GUID
#'
#' Get the data for any given GUID and return it as a data frame. Note that we
#' use the "ajson" json format from the API. The "json" format has a more
#' complex structure.
#'
#' We do use the json response to get the \code{fLength} value, which indicates the
#' length of the dataset. This way we can include a fixed way to get around the
#' default limit of 1000 rows of the Junar API.
#'
#' Note that this removes all meta-data from the json response given by the API.
#'
#' @param base_url The base URL of the Junar service
#' @param api_key The user's API key for the Junar service
#' @param guid The GUID of the data set of interest
#' @keywords GUID
#' @export

get_data <- function(base_url, api_key, guid) {
  if (missing(base_url)) {
    warning("Please add a valid base URL")
  }
  if (missing(api_key)) {
    warning("Please add a valid API key for the base URL you are trying to access")
  }
  if (missing(guid)) {
    warning("Please add a valid GUID for the dataset you are trying to access")
  }
  try({
    r_json <- GET(paste(base_url, guid, "/data.json/","?auth_key=", api_key, sep=""), accept_json())
    jsondata <- fromJSON(content(r_json, "text"))
    data_length <- jsondata$result$fLength

    r_ajson <- GET(paste(base_url, guid, "/data.ajson/","?auth_key=", api_key, "&limit=", data_length, sep=""), accept_json())
    dataset <- fromJSON(content(r_ajson, "text"))
    dataset <- dataset$result
    df <- as.data.frame(dataset)
    colnames(df) <- dataset[1,]
    df <- df [-1,]
    return(df)
  })
}
