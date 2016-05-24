#' Reset a monitor
#'
#' \code{uptimerobot.monitor.reset} remove all the statistics and logs associated to a monitor ID.
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
#' # Get a list of all available monitors, and take the first id
#' monitors.id <- uptimerobot.monitors(api.key, fields="id")[1,1]
#'  
#'  # Reset the stats for that monitor
#'  uptimerobot.monitor.reset(api.key, monitor.id)
#'  
#' }
#' @importFrom RCurl getURL
#' @importFrom rjson fromJSON
#' @importFrom utils URLencode
#' @export 
uptimerobot.monitor.reset <- function(api.key, id){
  
  if(is.null(api.key) | 
     is.na(api.key) | 
     (is.character(api.key) & nchar(api.key)==0)
  ) stop("api.key cannot be empty or NULL")
  
  data <- fromJSON(
    getURL(
      URLencode(paste0("https://api.uptimerobot.com/resetMonitor?apiKey=",
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
