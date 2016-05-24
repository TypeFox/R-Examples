#' @title GetFitDatasource
#' @rdname GetFitDatasource
#' @export
#'
#' @param token - OAuth 2.0 access token
#' @param datasource - Data Stream ID
#' @description
#' Returns a list with datasource attributes.
#' Refer to \url{https://developers.google.com/fit/rest/v1/data-types} for full documentation.

GetFitDatasource <- function(token,datasource) {

  return (fromJSON(getURL(URLencode(paste("https://www.googleapis.com/fitness/v1/users/me/dataSources/",
                                datasource,
                                sep="")
                                ),
                          httpheader = FitHTTPHeader(token))))

}
