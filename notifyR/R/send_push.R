send_push <-
function(user,message,title="",device="",url="",url_title="",priority=FALSE,timestamp=Sys.time()) {
  token <- "fY1h9Ph686pOtjUEsOk5hEZpkimzUu"
  if (!missing(user)) {
    } else {
      stop("User ID is missing.",call.=FALSE)
    }
  if (!missing(message)) {
    } else {
    stop("Message is missing.",call.=FALSE)
    }
  if (url.exists("https://api.pushover.net/",.opts=list(ssl.verifypeer = FALSE))) {
    } else {
    stop("No Connection to API.",call.=FALSE)
    }
 
  if (nchar(message)>512) {
    warning("The message has more than 512 characters and will be shortened.",call.=FALSE)
    message <- substr(message,1,512)
  }
  if (nchar(url)>500) {
    warning("The url has more than 500 characters and will be shortened.",call.=FALSE)
    message <- substr(url,1,500)
  }
  if (nchar(url_title)>50) {
    warning("The url title has more than 50 characters and will be shortened.",call.=FALSE)
    message <- substr(url_title,1,50)
  }
 if (!priority) {
   priority <- ""
   } else {
   priority <- "1"
   }
  
  po <- list(token = as.character(token), user = as.character(user), message = as.character(message), title = as.character(title), device = as.character(device), url = as.character(url), url_title = as.character(url_title), priority = as.character(priority), timestamp = as.character(floor(as.numeric(timestamp))))
  
  result <- suppressWarnings(postForm("https://api.pushover.net/1/messages.json",
                  .params=list(token = po$token,
                  user = po$user,
                  message = po$message,
                  title = po$title,
                  device = po$device,
                  url = po$url,
                  url_title = po$url_title,
                  priority = po$priority,
                  timestamp = po$timestamp),
                  .opts=list(ssl.verifypeer = FALSE)
                  ))
  result <- fromJSON(result)
  if (result$status == 1) {
    
  } else {
    variables <- character(length(result)-1)
    for (i in 1:length(variables)) {
      variables[i] <- names(result[i])
      }
    stop(sprintf(ngettext(length(variables),
                         "parameter %s is invalid\n",
                         "parameters %s are invalid\n"),
                paste(sQuote(variables), collapse=", ")),call.=FALSE)
  }
}
