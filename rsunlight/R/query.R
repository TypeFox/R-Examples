query <- function(url, args, ...){

  #Query
  query_results <- httr::GET(url, query = args, ...)

  #Check status. 403s are the most likely and we don't want to be just providing a shouty "403 FORBIDDEN";
  #we know why those show up and can be reasonable about it.
  status <- query_results$status_code
  if (status == 403) {
    stop("Your connection was not authorised: you have either not provided an API key, or provided
         an invalid one. Please see ?rsunlight", .call = FALSE)
  }

  #Check status against remaining options, check returned type.
  httr::stop_for_status(status)
  stopifnot(query_results$headers$`content-type` == 'application/json; charset=utf-8')

  #Assuming it hasn't blown up by now, return
  return(query_results)
}
