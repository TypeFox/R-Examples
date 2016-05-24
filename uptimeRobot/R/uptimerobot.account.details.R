#' Get the account details for who is linked to the given API key
#'
#' \code{uptimerobot.account.details} returns a list or a vector with the account details connected to the given api key.
#' 
#' @author Gabriele Baldassarre
#' 
#' @return A list or a vector with the account details.
#'
#' @param api.key string with a valid key for connecting to Uptimerobot API.
#' @param unlist logical. Set to \code{TRUE} to unlist the output to a named vector, \code{FALSE} to get a named list.
#'
#' @examples
#' \dontrun{
#' # Let's assume the api.key is available into the environment variable KEY
#' api.key <- Sys.getenv("KEY", "")
#' 
#' # Returns details as a list
#' details.list <- uptimerobot.account.details(api.key)
#' 
#' # Returns details as a vector
#' details.num <- uptimerobot.account.details(api.key, unlist = TRUE)
#' }
#' 
#' @importFrom RCurl getURL
#' @importFrom rjson fromJSON
#' @importFrom utils URLencode
#' @export
uptimerobot.account.details <- function(api.key, unlist = FALSE){
  
  if(is.null(api.key) | 
     is.na(api.key) | 
     (is.character(api.key) & nchar(api.key)==0)
  ) stop("api.key cannot be empty or NULL")
  
  data <- fromJSON(
    getURL(
      URLencode(paste0("https://api.uptimerobot.com/getAccountDetails?apiKey=",
                       api.key,
                       "&format=json&noJsonCallback=1"
      )      
      )
    ),
    unexpected.escape="keep"
  )
  
  if(data$stat=="ok") {
    
    if(!unlist) return(lapply(data$account, function(x){ as.integer(x)}))
    
    data.unlisted <- as.integer(unlist(data$account))
    names(data.unlisted) <- names(unlist(data$account))
    return(data.unlisted)
  }
  else {
    stop(data$message)
  }
  
}
