#' @title GetFitDataset
#' @rdname GetFitDataset
#' @export
#' @param token - OAuth 2.0 access token
#' @param datasource - Data Stream ID
#' @param startTime - Start time for dataset. \code{Date} datatype is required
#' @param endTime - End time for dataset. \code{Date} datatype is required
#' @description
#' Retrieves a dataset for a given \code{datasource} and time range defined by
#' \code{startTime} and \code{endTime}.
#' Refer to \url{https://developers.google.com/fit/rest/v1/datasets} for full documentation.

GetFitDataset <- function(token, datasource, startTime, endTime) {

  ValidateDateRange(startTime,endTime)
  ValidateDatasource(datasource)

  fitStartTS <- EpochTime(startTime)
  fitEndTS <- EpochTime(endTime)

  url <-
    paste(
      "https://www.googleapis.com/fitness/v1/users/me/dataSources/",
      datasource,
      "/datasets/",
      fitStartTS,
      "-",
      fitEndTS,
      sep = ""
    )

  url <- URLencode(url)
  
  print(paste("Requested URL:",url))
  
  return (fromJSON(getURL(url,
                          httpheader = FitHTTPHeader(token))))

}
