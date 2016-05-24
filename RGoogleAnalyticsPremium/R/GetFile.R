#' Query the Google Analytics Premium API for the specified dimensions, metrics and other query parameters
#'
#' @export
#'
#' @param query.builder Name of the object created using \code{\link{QueryBuilder}}
#'
#' @param token Name of the token object created using \code{\link{Auth}}
#'
#' @param accountid Google analytics premium account id
#' @param webpropertyid Webproperty ID for google analytics premium account
#' @param profileid View ID for google analytics premium account
#'
#' @return
#'    It returns path to file on local drive that contains extracted unsampled data.
#'
#' @importFrom httr status_code
#' @importFrom httr GET
#' @importFrom httr content
#' @importFrom jsonlite fromJSON
#' @importFrom jsonlite flatten
#'
GetFile <- function(query.builder, token, accountid, webpropertyid, profileid) {

  #' Code to fire query to the API
  query_url <- "https://www.googleapis.com/analytics/v3/management"
  url <- paste("accounts",accountid,"webproperties",webpropertyid,"profiles",profileid,"unsampledReports",sep = "/")
  query.uri <- paste(query_url,url,sep = "/")
  dataframe.param <- data.frame()

  # Set the CURL options for Windows
  options(RCurlOptions = list(capath = system.file("CurlSSL",
                                                   "cacert.pem",
                                                   package = "RCurl"),
                              ssl.verifypeer = FALSE))

  # Set all the Query Parameters
  query.builder$SetQueryParams()
  postbody <- NULL
  postbody <- ToBody(query.builder,token)

  res <- GetPostResponse(query.uri,postbody,token)

  if (status_code(res) == 200) {
     cat("Request to API successfully fired.\n")
  }
  else {
    cat("Request unsuccessfull.. please try again.\n")

  }

  #' Code to get the drive document id
  repeat
  {
    id <- content(res)$id
      if(!is.null(id))
    { break }
  }

  uri <- paste(query.uri,id,sep = "/")
  cat("Downloading of data is in progress.. ")
  repeat
  {
      load("token_files")
      ValidateToken(token)
      tok <- paste("Bearer", token$credentials$access_token)
      resp <- GET(uri,add_headers(Authorization = tok))
      data_json <- fromJSON(content(resp, as = "text"), simplifyVector = TRUE, flatten = flatten)
      if(!is.null(data_json$driveDownloadDetails$documentId))
        { break }
  }

  #print(data_json$driveDownloadDetails$documentId)

  #' Code to get the download url for drive document id
  doc_id <- data_json$driveDownloadDetails$documentId
  get_doc_uri <- "https://www.googleapis.com/drive/v2/files"
  uri <- paste(get_doc_uri,doc_id,sep="/")
  response <- GET(uri,add_headers(Authorization = tok))
  data_list <- fromJSON(content(response, as = "text"), simplifyVector = TRUE, flatten = flatten)


  #' Code to download data to the local drive
  url <- data_list$downloadUrl
  response_content <- GET(url,add_headers(Authorization = tok))
  data.frame <- content(response_content,"text")
  file_nm <- data_json$title
  ext <- "csv"
  path <- paste(file_nm,ext,sep = ".")
  write(data.frame,file = path,append = TRUE, sep = "\n")
  cat("Your data is downloaded.")
  return(path)
}
