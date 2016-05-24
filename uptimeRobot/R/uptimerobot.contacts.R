#' Get general informations about the alert contacts
#'
#' \code{uptimerobot.contacts} extracts a dataset with general informations
#' for a set of contacts used to be alert in case of up/down of the given monitors.
#' 
#' @details 
#' The alert contacts are whom to be notified when the monitor goes up/down. 
#' 
#' If a vector of contact IDs is not given, the function will return data for all the available contacts.
#' 
#' The API uses pagination and returns no more than 50 contacts on each page. Use \code{limit} and \code{offset} to set a different number of
#' monitors to get on each page and to move between pages. Leave default values to get all the data.
#' 
#' @return A dataset with general informations about the contacts.
#' 
#' @author Gabriele Baldassarre
#' 
#' @seealso \code{\link{uptimerobot.monitors}}
#'
#' @param api.key A valid key for connecting to UptimeRobors public API.
#' @param contacts vector or comma-delimited string with the IDs of the contacts to get.
#' If the argument is NULL or missing, all the available contacts will be returned. 
#' @param fields vector or comma-delimited string with the general informations to include in the output dataset.
#' You may use the helper function \code{\link{uptimerobot.fields}} if you don't want to manually compile the list of fields.
#' @param limit An integer value used for pagination. Defines the max number of records to return in each page. Default and max. is 50.
#' @param offset An integer value to set the index of the first monitor to get (used for pagination).
#'
#' @examples
#' \dontrun{
#' # Let's assume the api.key is available into the environment variable KEY
#' api.key <- Sys.getenv("KEY", "")
#' 
#' # Returns all the contacts with a default set of attributes
#' contacts.df <- uptimerobot.contacts(api.key)
#' 
#' # Returns all the contacts and all the attributes
#' contacts.full.df <- uptimerobot.contacts(api.key, fields=uptimerobot.fields("contact")$full))
#' 
#' # Returns only the two contacts with ID: 1234, 5678
#' contacts.df <- uptimerobot.contacts(api.key, c("1234", "5678"))
#' }
#' 
#' @importFrom RCurl getURL
#' @importFrom rjson fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom utils URLencode
#' @export 
uptimerobot.contacts <- function(api.key, 
                                  contacts=NULL,
                                  limit=50,
                                  offset=0,
                                  fields=uptimerobot.fields("contact")$typical){
  
  fields <- (function(x){ if(length(x) == 0) "id" else x })(fields)

  fields.o <- unique(unlist(strsplit(fields, split = ",")))
  fields.v <- unique(c(unlist(strsplit(fields, split = ","))), "id")

  data <- fromJSON(
    getURL(
      URLencode(paste0("https://api.uptimerobot.com/getAlertContacts?apiKey=",
             api.key,
             ifelse(is.null(contacts), "", paste0("&alertcontacts=", paste0(unique(unlist(strsplit(contacts, split = ","))), collapse = "-"), sep="")),
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
        rbind.fill,lapply(data$alertcontacts$alertcontact, function(x) {
          return(do.call(data.frame, list(x[which(names(x) %in%  fields.v)], stringsAsFactors = FALSE)))
        })
      )
      
      # Convert to proper datatypes
      if("type" %in% fields.o) data.merged$type <- factor(as.integer(data.merged$type), levels=c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11), labels=c("SMS", "Email", "Twitter DM", "Boxcar", "WebHook", "Pushbullet", "Zapier", "Pushover", "Hipchat", "Slack"))
      if("status" %in% fields.o) data.merged$status <- factor(as.integer(data.merged$status), levels=c(0, 1, 2), labels=c("not activated", "paused", "active"))
      
      if(!("id" %in% fields.o)) data.merged$id <- NULL
      
      # Pagination
      if((as.integer(data$offset) + as.integer(data$limit)) >= as.integer(data$total)) return(data.merged)
      else {
        rbind.fill(data.merged, uptimerobot.contacts(api.key = api.key, 
                                                     contacts = contacts, 
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
