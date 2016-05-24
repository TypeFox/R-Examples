getAllTicketMetrics <- function(){
  curl <- getCurlHandle()
  result <- list()
  stopPaging <- FALSE
  i <- 1
  
  ## Need to page through the results since only 100 are returned at a time
  while(stopPaging==FALSE){
    result[[i]]<-getURL(paste(.ZendeskEnv$data$url, .ZendeskEnv$data$ticket_metrics, "?page=", i, sep=""), curl=curl, ssl.verifypeer=FALSE,
                        .opts=list(userpwd=(paste(.ZendeskEnv$data$username, .ZendeskEnv$data$password, sep=":"))))
    if(is.null(fromJSON(result[[i]])$next_page)){
      stopPaging = TRUE
    }
    i <- i + 1
  }
  
  ## Transform the JSON data to a data.frame
  json.data <- lapply(unlist(result), fromJSON)
  pre.result <- lapply(json.data, function(x) do.call("rbind", x$ticket_metrics))
  final.result<-do.call("rbind", pre.result)
  ticket_metrics.df <- data.frame(final.result)
  ticket_metrics.df <- unlistDataFrame(ticket_metrics.df)
  return(ticket_metrics.df)
}        

