#'@importFrom httr GET stop_for_status content user_agent
threeword_query <- function(url, ...){
  result <- httr::GET(url, httr::user_agent("what3words R client - https://www.github.com/ironholds/threewords/"), ...)
  httr::stop_for_status(result)
  output <- httr::content(x = result, as = "parsed")
  if("error" %in% names(output)){
    stop(output$message)
  }
  return(output)
}

clean <- function(result){
  result$words <- unlist(result$words)
  result$position <- as.numeric(unlist(result$position))
  return(result)
}