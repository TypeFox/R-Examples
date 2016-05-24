#'@importFrom httr stop_for_status GET user_agent content
pv_query <- function(params, ...){
  url <- paste0("https://wikimedia.org/api/rest_v1/metrics/pageviews/", params)
  result <- httr::GET(url, httr::user_agent("pageviews API client library - https://github.com/Ironholds/pageviews"))
  
  # Specific check for missing data
  if(result$status_code == 404){
    stop("The pageview data available does not cover the range you requested")
  }
  
  # General check and return
  httr::stop_for_status(result)
  return(httr::content(result))
}
