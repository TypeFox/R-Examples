#'@title Create a token for WhoAPI
#'@description This function allows the creation of connectors to whoapi, which can be passed
#'between function calls, providing consistent API keys and options for different calls.
#'
#'@param key the API key to use.
#'
#'@param user_agent the user agent to use. NULL by default, which will use the \code{whoapi} package user agent.
#'
#'@export
whoapi_token <- function(key, user_agent = NULL){

  if(is.null(user_agent)){
    user_agent <- "whoapi R client: https://github.com/ironholds/whoapi"
  }

  output_object <- list(key = key, user_agent = user_agent)

  return(output_object)
}


#Convert character-stored numeric values to logical
char_convert <- function(x){
  return(as.logical(as.numeric(x)))
}

#Query function
whoapi_query <- function(token, url, ...){
  url <- paste0("http://api.whoapi.com/?apikey=", token$key, url)
  result <- httr::GET(url, user_agent(token$user_agent), ...)
  httr::stop_for_status(result)
  result <- httr::content(result)
  if(char_convert(result$status)){
    stop(result$status_desc)
  }
  return(result)
}
