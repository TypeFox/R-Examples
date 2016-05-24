#' Get contacts informations for one or more monitors
#'
#' \code{uptimerobot.monitor.contacts} return a dataset with all the alert contacts that will be triggered
#' when a log is collected for the given monitors IDs.
#' 
#' The API uses pagination and returns no more than 50 monitors on each page. Use \code{limit} and \code{offset} to set a different number of
#' monitors to get on each page and to move between pages. Leave default values to get all the data.
#' 
#' @return A dataset with the alert contacts.
#' 
#' @author Gabriele Baldassarre
#' 
#' @seealso \code{\link{uptimerobot.monitors}}, \code{\link{uptimerobot.monitor.logs}}, \code{\link{uptimerobot.monitor.responses}}
#' 
#' @param api.key A valid key for connecting to UptimeRobors public API.
#' @param monitors vector or comma-delimited string with the IDs of the monitors to get.
#' @param limit An integer value used for pagination. Defines the max number of records to return in each page. Default and max. is 50.
#' @param offset An integer value to set the index of the first monitor to get (used for pagination).
#' 
#' @examples
#' \dontrun{
#' # Let's assume the api.key is available into the environment variable KEY
#' api.key <- Sys.getenv("KEY", "")
#' 
#' # Returns all the monitors IDs. Since the function always return a data.frame
#' # (even if you ask only for a column), you have to reference the column to get a character vector.
#' monitors.id <- uptimerobot.monitors(api.key, fields="id")$id
#' 
#' # Returns all the contacts registered for the given monitors
#' logs.df <- uptimerobot.monitor.contacts(api.key, monitors=monitors.id)
#' }
#' 
#' @importFrom RCurl getURL
#' @importFrom rjson fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom utils URLencode
#' @export 
uptimerobot.monitor.contacts <- function(api.key, 
                                         monitors,
                                         limit=50,
                                         offset=0){
  
  data <- fromJSON(
    getURL(
      URLencode(paste0("https://api.uptimerobot.com/getMonitors?apiKey=",
             api.key,
             "&monitors=", paste0(unique(unlist(strsplit(monitors, split = ","))), collapse = "-"),
             "&showMonitorAlertContacts=1",
             "&limit=", limit,
             "&offset=", offset,
             "&format=json&noJsonCallback=1"
      )      
    )
    ),
    unexpected.escape="keep"
  )
  
  if(data$stat=="ok") {
    return((function() {
      data.merged <- do.call(
        rbind.fill,lapply(data$monitors$monitor, function(x) {
          
          contacts <- do.call(rbind, lapply(x$alertcontact, function(y){
            do.call(data.frame, list(y, stringsAsFactors = FALSE))
          }))
          contacts$monitor.id <- x$id
          
          # Convert to proper datatypes
          contacts$type <- factor(as.integer(contacts$type), levels=c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11), labels=c("SMS", "Email", "Twitter DM", "Boxcar", "WebHook", "Pushbullet", "Zapier", "Pushover", "Hipchat", "Slack"))
          
          
          return(contacts[,c(6,1,2,3,4,5)])
        })
      )
      
      # Convert to proper datatypes
      data.merged$threshold <- as.integer(data.merged$threshold)
      data.merged$recurrence <- as.integer(data.merged$recurrence)
      
      # Pagination
      if((as.integer(data$offset) + as.integer(data$limit)) >= as.integer(data$total)) return(data.merged)
      else {
        rbind.fill(data.merged, uptimerobot.monitor.contacts(api.key = api.key, 
                                                             monitors = monitors, 
                                                             limit = limit,
                                                             offset = as.character(as.integer(offset) + as.integer(limit))))
      }
    })()
    )
  }
  else {
    stop(data$message)
  }
  
}
