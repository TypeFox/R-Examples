user.info <-
function(token=NA,consumer_key=NA,consumer_secret=NA){
  
  if(class(token)[1]!="Token1.0")
    stop("token must be a Token1.0 type")
  
  if(!is.character(consumer_key))
    stop("consumer_key must be a string")
  
  if(!is.character(consumer_secret))
    stop("consumer_secret must be a string")
  
  url<-"http://api.tumblr.com/v2/user/info"
  connection="GET"
  bodyParams<-list()
  res<-fromJSON(http.connection(url,token,bodyParams,consumer_key,consumer_secret,connection))
  
  return(res)
}
