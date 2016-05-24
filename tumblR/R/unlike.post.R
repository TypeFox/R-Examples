unlike.post <-
function(id=NA,reblog_key=NA,token=NA,consumer_key=NA,consumer_secret=NA){
  
  if(class(token)[1]!="Token1.0")
    stop("token must be a Token1.0 type")
  
  if(!is.character(consumer_key))
    stop("consumer_key must be a string")
  
  if(!is.character(consumer_secret))
    stop("consumer_secret must be a string")
  
  if(is.na(id)){
    stop("id is a required field")
  } else{
    if(!is.numeric(id))
      stop("id must be a numeric type")
  }
  
  if(is.na(reblog_key)){
    stop("reblog_key is a required field")
  } else{
    if(!is.character(reblog_key))
      stop("reblog_key must be a string type")
  }
  
  url<-"http://api.tumblr.com/v2/user/unlike"
  bodyParams <- list(id=as.character(id), reblog_key=reblog_key)
  connection<-"POST"
  
  res<-fromJSON(http.connection(url,token,bodyParams,consumer_key,consumer_secret,connection))
  
  return(res)
}
