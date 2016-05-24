post.delete <-
function(base_hostname=NA,id=NA,token=NA,consumer_key=NA,consumer_secret=NA){

  if(!is.character(base_hostname))
    stop("base_hostname must be a string")
  
  if(class(token)[1]!="Token1.0")
    stop("token must be a Token1.0 type")
  
  if(!is.character(consumer_key))
    stop("consumer_key must be a string")
  
  if(!is.character(consumer_secret))
    stop("consumer_secret must be a string")
  
  if(is.na(id)){
    stop("id is a requested field")
    
  } else {
    if(!is.numeric(id))
      stop("id must be a number")
  }
  
  url<-paste("http://api.tumblr.com/v2/blog/",base_hostname,"/post/delete",sep="")
  bodyParams <- list(id=as.character(id))
  connection<-"POST"
  
  res<-fromJSON(http.connection(url,token,bodyParams,consumer_key,consumer_secret,connection))
  
  return(res)
}
