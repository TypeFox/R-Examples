unfollow <-
function(url=NA,token=NA,consumer_key=NA,consumer_secret=NA){
  
  if(is.na(url)){
    stop("url is a required field")
  } else{
    if(!is.character(url))
      stop("url must be a string type")
  }
  
  if(class(token)[1]!="Token1.0")
    stop("token must be a Token1.0 type")
  
  if(!is.character(consumer_key))
    stop("consumer_key must be a string")
  
  if(!is.character(consumer_secret))
    stop("consumer_secret must be a string")
  
  url_post<-"http://api.tumblr.com/v2/user/unfollow"
  bodyParams <- list(url=url)
  connection<-"POST"
  
  res<-fromJSON(http.connection(url_post,token,bodyParams,consumer_key,consumer_secret,connection))
  
  return(res)
}
