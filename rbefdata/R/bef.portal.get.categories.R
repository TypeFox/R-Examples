#' Get a data groups categories as data frame (data admin interface)
#'
#' This function gets you a data groups categories as data frame. You can use the
#' data frame to harmonize your data using the last column and the IDs in the
#' first row. After you have relinked the data you can upload the data frame using
#' the categories upload command.
#'
#' @param datagroup_id It the ID of the data group you like to get the categories from.
#' @param curl Pass in a curl handle here if you like with own options or to reduce the
#'        memory footprint.
#' @param ... Additional parameters passed to curl.
#' @return Returns a data frame with the categories ID, short and long description and
#          merge ID.
#' @examples \dontrun{
#'   bef.portal.get.categories_for(datagroup_id = 22)
#'       }
#' @import RCurl
#' @export bef.portal.get.categories_for bef.get.categories_for
#' @aliases bef.get.categories_for

bef.portal.get.categories_for <- bef.get.categories_for <- function(datagroup_id, curl = getCurlHandle(), ...) {
  datagroup_url = datagroups_url(datagroups_id = datagroup_id, type = "download")
  response_body = getURLContent(datagroup_url, curl = curl, ...)
  if(getCurlInfo(curl)$response.code != 200) {
    msg = sprintf("Datagroups (id=%d) not found or not accessible. Please check your credentials!", id)
    stop(msg)
  }
  dataframe = read.csv(text = response_body)
  return(dataframe)
}
