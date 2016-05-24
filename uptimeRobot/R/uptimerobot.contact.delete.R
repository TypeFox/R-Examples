#' Delete an alert contact
#'
#' \code{uptimerobot.contact.delete} removes an alert contanct, unlinking from all the registered monitors.
#' 
#' @return The function returns \code{TRUE} in case success. An error is thrown otherwise.
#'  
#' @author Gabriele Baldassarre
#' 
#' @param api.key string with a valid key for connecting to Uptimerobot API.
#' @param id numeric or integer with the ID of the contact to delete.
#'
#' @examples
#' \dontrun{
#'  # Let's assume the api.key is available into the environment variable KEY
#'  api.key <- Sys.getenv("KEY", "")
#'  
#'  # Delete the contact with id=12345678
#'  if(uptimerobot.contact.delete(api.key, 12345678){
#'    message("Alert contact successfully deleted!")
#'  }
#' }
#' 
#' @importFrom RCurl getURL
#' @importFrom rjson fromJSON
#' @importFrom utils URLencode
#' @export 
uptimerobot.contact.delete <- function(api.key, id){
  
  if(is.null(api.key) | 
     is.na(api.key) | 
     (is.character(api.key) & nchar(api.key)==0)
  ) stop("api.key cannot be empty or NULL")
  
  data <- fromJSON(
    getURL(
      URLencode(paste0("https://api.uptimerobot.com/deleteAlertContact?apiKey=",
             api.key,
             "&alertcontactID=", id,
             "&format=json&noJsonCallback=1"
      )      
    )
    ),
    unexpected.escape="keep"
  )
  
  if(data$stat=="ok") {
    return(TRUE)
  }
  else {
    stop(data$message)
  }
  
}
