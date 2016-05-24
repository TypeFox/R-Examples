#' @title GetFitDatasources
#' @rdname GetFitDatasources
#' @export
#' @description
#' \code{GetFitDatasources} returns a dataframe with user's datasources.
#' @param token - OAuth 2.0 access token

GetFitDatasources <- function(token) {

  return (fromJSON(getURL("https://www.googleapis.com/fitness/v1/users/me/dataSources",
                 httpheader = FitHTTPHeader(token)))[[1]])

}
