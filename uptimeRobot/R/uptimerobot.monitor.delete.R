#' Delete a monitor
#'
#' \code{uptimerobot.monitor.delete} remove a monitor and all existing statistics of it.
#' 
#' @return
#' The function returns \code{TRUE} in case success. An error is thrown otherwise.
#'  
#' @author Gabriele Baldassarre
#' 
#' @param api.key string with a valid key for connecting to Uptimerobot API.
#' @param id numeric or integer with the ID of the monitor to delete.
#'
#' @examples
#' \dontrun{
#'  # Let's assume the api.key is available into the environment variable KEY
#'  api.key <- Sys.getenv("KEY", "")
#'  
#'  # Create a monitor and get its monitor.id
#'  monitor.id <- uptimerobot.monitor.new(api.key, 
#'    friendly.name="Open Analytics", 
#'    url="https://gabrielebaldassarre.com", type="http"
#'  )
#'  
#'  # Change the friendly name of the monitor
#'   if(uptimerobot.monitor.edit(api.key, 
#'      monitor.id, 
#'      friendly.name="Open Analytics - gabrielebaldassarre.com"
#'   ){
#'    message("Monitor has been successfully edited!")
#'  }
#'  
#'  # Delete the just-made monitor
#'  if(uptimerobot.monitor.delete(api.key, monitor.id){
#'    message("Monitor has been successfully deleted!")
#'  }
#' }
#' 
#' @importFrom RCurl getURL
#' @importFrom rjson fromJSON
#' @importFrom utils URLencode
#' @export 
uptimerobot.monitor.delete <- function(api.key, id){
  
  if(is.null(api.key) | 
     is.na(api.key) | 
     (is.character(api.key) & nchar(api.key)==0)
  ) stop("api.key cannot be empty or NULL")
  
  data <- fromJSON(
    getURL(
      URLencode(paste0("https://api.uptimerobot.com/deleteMonitor?apiKey=",
                       api.key,
                       "&monitorID=", id,
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
