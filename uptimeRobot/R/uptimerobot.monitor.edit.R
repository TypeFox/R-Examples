#' Edit a monitor
#'
#' \code{uptimerobot.monitor.edit} edits the properties for an existing monitor.
#' 
#' @details
#' If a property has not to be updated, just omit it from the parameters or set to \code{NA}.
#' To erase the value of a property, set it to an empty string, ie \code{""}, instead (not \code{NA} or \code{NULL}!).
#' 
#' The type of a monitor can not be edited (like changing a HTTP monitor into a Port monitor).
#'  
#' The alert contacts are whom to be notified when the monitor goes up/down.
#' 
#' Multiple alert contact IDs can be sent in a character vector or in a data frame. If you pass alert contact IDs in a vector, each element must be formatted in the form \code{<id>_<threshold>_<recurrence>} (note the underscores).
#' If you prefer to format it as a data.frame, it must have these three columns: \code{id, threshold, recurrence}, numeric or integer. Order of the columns doesn't matter.
#' 
#' Please note that thresholds and recurrences can be omitted (default to zero) and, as they are only available in the Pro Plan, they are always 0 in the Free Plan.
#'
#' @return 
#' The function returns \code{TRUE} in case success. An error is thrown otherwise.  
#'  
#' @author Gabriele Baldassarre
#' 
#' @param api.key string with a valid key for connecting to Uptimerobot API.
#' @param id numeric or integer with the ID of the monitor to edit.
#' @param friendly.name string the friendly (screen) name of the monitor.
#' @param url string with the URL/IP of the monitor.
#' @param activate logical to set the status of the monitor. Set to \code{TRUE} to start the monitor or \code{FALSE} to put it in paused state.
#' @param subtype string used only for "Port monitoring" to set which pre-defined port/service is monitored or if a custom port is monitored. You can use both the friendly name (string) or the index (integer) here.
#' @param port string used only for "Port monitoring" to set the port monitored.
#' @param keyword.type required string in Keyword monitoring".
#' @param keyword.value string with the value of the keyword (required for keyword monitoring).
#' @param HTTP.username string used for password-protected web pages (HTTP Basic Auth). Set to empty string to erase the current username. Available for HTTP and keyword monitoring.
#' @param HTTP.password string used for password-protected web pages (HTTP Basic Auth). Set to empty string to erase the current password.Available for HTTP and keyword monitoring.
#' @param alert.contacts character vector or data frame with the IDs to alert each with their threshold and recurrence values.
#' @param interval integer with the interval for the monitoring check (in minutes).
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
uptimerobot.monitor.edit <- function(api.key,
                                     id,
                                     friendly.name=NULL,
                                     url=NULL,
                                     activate=TRUE,
                                     subtype=NULL,
                                     port=NULL,
                                     interval=NULL,
                                     keyword.type=NULL,
                                     keyword.value=NULL,
                                     HTTP.username=NULL,
                                     HTTP.password=NULL,
                                     alert.contacts=NULL){
  
  if(is.null(api.key) | 
     is.na(api.key) | 
     (is.character(api.key) & nchar(api.key)==0)
  ) stop("api.key cannot be empty or NULL")
  
  # Decode monitor subtype
  if(!(is.null(subtype))){
    if(class(subtype) == "character"){
      subtype <- as.numeric(factor(toupper(subtype), labels=c(1,2,3,4,5,6,99), levels=c("HTTP", "HTTPS", "FTP", "SMTP", "POP3", "IMPAP", "Custom Port")))
    } else if(!(class(subtype) %in% c("integer", "numeric"))) stop(paste0(class(subtype), "is not a valid format for monitor subtype", sep=" "))
  }
  
  if(!(is.null(alert.contacts))){
    if(is.data.frame(alert.contacts)) {
      
      if(!("threshold" %in% names(alert.contacts))) alert.contacts$threshold <- 0
      if(!("recurrence" %in% names(alert.contacts))) alert.contacts$recurrence <- 0
      
      alert.contacts <- paste(alert.contacts$id, alert.contacts$threshold, alert.contacts$recurrence, sep="_")
    }
    alert.contacts <- paste0(paste(alert.contacts$id, ifelse(is.na(alert.contacts$threshold), 0, alert.contacts$threshold), ifelse(is.na(alert.contacts$recurrence), 0, alert.contacts$recurrence), sep="_"), collapse="-")
    
  }
  
  data <- fromJSON(
    getURL(
      URLencode(paste0("https://api.uptimerobot.com/editMonitor?apiKey=",
                       api.key,
                       "&monitorID=", id,
                       ifelse(is.null(friendly.name), "", paste0("&monitorFriendlyName", friendly.name, sep="=")),
                       ifelse(is.null(url), "", paste0("&monitorURL", url, sep="=")),
                       "&monitorStatus=", as.integer(activate),
                       ifelse(is.null(subtype), "", paste0("&monitorSubType", subtype, sep="=")),
                       ifelse(is.null(port), "", paste0("&monitorPort", port, sep="=")),
                       ifelse(is.null(interval), "", paste0("&monitorInterval", interval, sep="=")),
                       ifelse(is.null(keyword.type), "", paste0("&monitorKeywordType", keyword.type, sep="=")),
                       ifelse(is.null(keyword.value), "", paste0("&monitorKeywordValue", keyword.value, sep="=")),
                       ifelse(is.null(HTTP.username), "", paste0("&monitorHTTPUsername", HTTP.username, sep="=")),
                       ifelse(is.null(HTTP.password), "", paste0("&monitorHTTPPassword", HTTP.password, sep="=")),
                       ifelse(is.null(alert.contacts), "", paste0("&monitorAlertContacts", alert.contacts, sep="=")),
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
