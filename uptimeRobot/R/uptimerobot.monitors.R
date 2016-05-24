#' Get general informations about monitors
#'
#' \code{uptimerobots.monitors.responses} return a dataset with general informations
#' for a set of monitors.
#' 
#' @details 
#' If a vector of monitor is not given, the function will return data for all the available monitors.
#' 
#' \code{summary} parameter expect a lists of three named logic values that set which columns of additional statistics for each monitor must  be added to output dataset for
#' each available monitor. These are summary values only, as the instances are obtained using a set of dedicated functions. 
#' 
#' \enumerate{
#'   \item \code{response.times} set to \code{TRUE} to add a column with the number of pings with response times available for the monitor to the output. These values can be queried using \code{\link{uptimerobot.monitor.responses}} function.
#'   \item \code{log.records} set to \code{TRUE} to add a column with the number of log entries recorded for the monitor to the output. These records can be queried using \code{\link{uptimerobot.monitor.logs}} function.
#'   \item \code{alert.contacts} set to \code{TRUE} to add a column with the number of alert contacts binded to the monitor to the output. Detailed informations about these contacts can be queried using \code{\link{uptimerobot.monitor.contacts}} function.
#' }
#' 
#' You may just add the elements you want to include into the list, as they default to \code{FALSE} if missing. Set an empty list to exclude all the summary statistics from the output.
#'
#' The API uses pagination and returns no more than 50 monitors on each page. Use \code{limit} and \code{offset} to set a different number of
#' monitors to get on each page and to move between pages. Leave default values to get all the data.
#' 
#' @return
#' A dataset with general informations about the given monitors
#' 
#' @author Gabriele Baldassarre
#' @seealso \code{\link{uptimerobot.monitor.responses}}, \code{\link{uptimerobot.monitor.logs}}, \code{\link{uptimerobot.monitor.contacts}}
#'
#' @param api.key A valid key for connecting to UptimeRobors public API.
#' @param monitors vector or comma-delimited string with the IDs of the monitors to get. If not used or set to \code{NULL}, will return all monitors in an account.
#' @param summary list of logical values to flag summary indicators to add to the output dataset.
#' @param types vector or comma-delimited string of monitor types. If not used or set to \code{NULL}, the function will return all monitors types (HTTP, keyword, ping..) in an account. Else, it is possible to define any number of monitor types. You can use both the friendly name (string) or the index (integer) here.
#' @param statuses vector or comma-delimited string of monitor statuses. If not used or set to \code{NULL}, the function will return  all monitors statuses (up, down, paused) in an account. Else, it is possible to define any number of monitor statuses. You can use both the friendly name (string) or the index (integer) here.
#' @param search An optional keyword of to search within monitor URL or friendly name to get filtered results.
#' @param limit An integer value used for pagination. Defines the max number of records to return in each page. Default and max. is 50.
#' @param offset An integer value to set the index of the first monitor to get (used for pagination).
#' @param fields vector or comma-delimited string with the general informations to include in the output dataset.
#' You may want to use the helper function \code{\link{uptimerobot.fields}} if you don't want to manually compile the list of fields.
#'
#' @examples
#' \dontrun{
#' # Let's assume the api.key is available into the environment variable KEY
#' api.key <- Sys.getenv("KEY", "")
#' 
#' # Returns all the  monitors with a default set of attributes
#' monitors.df <- uptimerobot.monitors(api.key)
#' 
#' #' # Returns all the monitors of 'keyword' type
#' monitors.kwd..df <- uptimerobot.monitors(api.key, type="keyword")
#' 
#' # Returns all the monitors and all the attributes
#' monitors.full.df <- uptimerobot.monitors(api.key, fields=uptimerobot.fields("monitor")$full))
#' 
#' # Returns only the two monitors with ID: 1234, 5678
#' monitors.df <- uptimerobot.monitors(api.key, c("1234", "5678"))
#' }
#' 
#' @importFrom RCurl getURL
#' @importFrom rjson fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom utils URLencode
#' @export 
uptimerobot.monitors <- function(api.key, 
                                 monitors=NULL,
                                 types=NULL,
                                 statuses=NULL,
                                 search=NULL,
                                 summary=list(),
                                 limit=50,
                                 offset=0,
                                 fields=uptimerobot.fields("monitor")$typical){
  
  if(is.null(api.key) | 
     is.na(api.key) | 
     (is.character(api.key) & nchar(api.key)==0)
  ) stop("api.key cannot be empty or NULL")
  
  fields <- (function(x){ if(length(x) == 0) "id" else x })(fields)
  
  fields.o <- unique(unlist(strsplit(fields, split = ",")))
  fields.v <- unique(c(unlist(strsplit(fields, split = ","))), "id")
  
  include.stat <- function(element){
    if(element %in% names(summary)){
      return(as.logical(summary[element]))
    } else return(FALSE)
  }
  include.responses <- include.stat("response.times")
  include.logs <- include.stat("log.records")
  include.contacts <- include.stat("alert.contacts")
  
  # Decode monitor types
  if(!(is.null(types))){
    if(class(types) == "character"){
      types <- unlist(strsplit(types, split = ","))
      if(suppressWarnings(all(is.na(as.numeric(types))))) types <- as.numeric(as.character(factor(tolower(types), labels=c("1", "2", "3", "4"), levels=c("http", "keyword", "ping", "port"))))
      else types <- as.numeric(types)
    } else if(!(class(types) %in% c("integer", "numeric"))) stop(paste0(class(types), "cannot be coerced to express a monitor status type", sep=" "))
  }
  
  # Decode monitor statuses
  if(!(is.null(statuses))){
    if(class(statuses) == "character"){
      statuses <- unlist(strsplit(statuses, split = ","))
      if(suppressWarnings(all(is.na(as.numeric(statuses))))) statuses <- as.numeric(as.character(factor(tolower(statuses), labels=c("0", "1", "2", "8", "9"), levels=c("paused", "not checked yet", "up", "seems down", "down"))))
      else statuses <- as.numeric(statuses)
    } else if(!(class(statuses) %in% c("integer", "numeric"))) stop(paste0(class(statuses), "cannot be coerced to express a monitor status attribute", sep=" "))
  }
   
  data <- fromJSON(
    getURL(
      URLencode(paste0("https://api.uptimerobot.com/getMonitors?apiKey=",
             api.key,
             ifelse(is.null(monitors), "", paste0("&monitors=", paste0(unique(unlist(strsplit(monitors, split = ","))), collapse = "-"), sep="")),
             ifelse(!is.null(types), paste0("&types=", paste0(unique(types[!is.na(types)]), collapse = "-"), sep=""), ""),
             ifelse(!is.null(statuses), paste0("&statuses=", paste0(unique(statuses[!is.na(statuses)]), collapse = "-"), sep=""), ""),
             ifelse(is.null(search), "", paste0("&search=", search, sep="")),
             "&responseTimes=", as.integer(include.responses),
             "&logs=", as.integer(include.logs),
             "&showMonitorAlertContacts=", as.integer(include.contacts),
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
          header <- do.call(data.frame, list(x[which(names(x) %in%  fields.v)], stringsAsFactors = FALSE))
          
          # Baking out the dataset
          if(include.responses) header$responses <- length(x$responsetime)
          if(include.logs) header$logs <- length(x$log)
          if(include.contacts) header$contacts <- length(x$alertcontact)
          
          return(header)
        })
      )
      
      # Convert to proper datatypes
      if("port" %in% fields.o) data.merged$port <- as.integer(data.merged$port)
      if("interval" %in% fields.o) data.merged$interval <- as.integer(data.merged$interval)
      
      # Lookup factors
      if("type" %in% fields.o) data.merged$type <- factor(as.integer(data.merged$type), levels=1:4, labels=c("HTTP", "Keyword", "Ping", "Port"))
      if("status" %in% fields.o) data.merged$status <- factor(as.integer(data.merged$status), levels=c(0, 1, 2, 8, 9), labels=c("paused", "not checked yet", "up", "seems down", "down"))
      if("subtype" %in% fields.o) data.merged$subtype <- factor(as.integer(data.merged$subtype), levels=c(1,2,3,4,5,6,99), labels=c("HTTP", "HTTPS", "FTP", "SMTP", "POP3", "IMPAP", "Custom Port"))
      if("keywordtype" %in% fields.o) data.merged$keywordtype <- factor(as.integer(data.merged$keywordtype), levels=c(1,2), labels=c("exists", "not exists"))
      
      
      if(!("id" %in% fields.o)) data.merged$id <- NULL
      
      # Pagination
      if((as.integer(data$offset) + as.integer(data$limit)) >= as.integer(data$total)) return(data.merged)
      else {
        rbind.fill(data.merged, uptimerobot.monitors(api.key = api.key, 
                                                     monitors = monitors,
                                                     types = types,
                                                     statuses = statuses,
                                                     summary = summary , 
                                                     limit = limit,
                                                     offset = as.character(as.integer(offset) + as.integer(limit)),
                                                     fields = fields)
        )
      }
    })()
    )
  }
  else {
    stop(data$message)
  }
  
}
