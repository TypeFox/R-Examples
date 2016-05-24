getAllSatisfactionRatings <- function(){
  curl = getCurlHandle()
  stopPaging <- FALSE
  result <- list()
  i <- 1
  
  ## Need to page through the results since only 100 are returned at a time
  while(stopPaging == FALSE){
    result[[i]]<-getURL(paste(.ZendeskEnv$data$url, .ZendeskEnv$data$satisfaction_ratings, "?page=" ,i, sep=""), curl=curl, ssl.verifypeer=FALSE,
                        .opts=list(userpwd=(paste(.ZendeskEnv$data$username, .ZendeskEnv$data$password, sep=":"))))    
    if(is.null(fromJSON(result[[i]])$next_page)){
      stopPaging <- TRUE
    }
    i <- i + 1
  }
  
  ## Transform the JSON data to a data.frame
  json.data <- lapply(unlist(result), fromJSON)
  pre.result <- lapply(json.data, function(x) do.call("rbind", x$satisfaction_ratings))
  final.result<-do.call("rbind", pre.result)
  satisfaction_ratings.df <- data.frame(final.result)
  satisfaction_ratings.df <- unlistDataFrame(satisfaction_ratings.df)
  return(satisfaction_ratings.df)
}
