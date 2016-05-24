#'@importFrom httr stop_for_status GET content user_agent
guardian_query <- function(path, ...){
  url <- paste0("https://content.guardianapis.com/", path)
  result <- httr::GET(url)
  httr::stop_for_status(result)
  content <- httr::content(result)
  if(length(content$response$results) == 0){
    stop("No results found")
  }
  return(content)
}

date_convert <- function(date){
  
  if(any(c("POSIXlt", "POSIXct") %in% class(date))){
    date <- as.Date(date)
  }
  return(as.character(date))
}

# Merges multiple args
merge_multis <- function(args){
  return(paste0(args, collapse = ","))
}
