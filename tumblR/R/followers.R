followers <-
function(base_hostname=NA,limit=20,offset=0,token=NA,consumer_key=NA,consumer_secret=NA){
  
  if(!is.character(base_hostname))
    stop("base_hostname must be a string")
  
  if(!is.numeric(limit) || (limit<=0 || limit>=21) )
    stop("limit must be a numeric type beetwen 0 and 20 (inclusive")
  
  if(!is.numeric(offset) || (offset<0 || offset>=21) )
    stop("offset must be a numeric type greater or equal to limit")
  
  if(class(token)[1]!="Token1.0")
    stop("token must be a Token1.0 type")
  
  if(!is.character(consumer_key))
    stop("consumer_key must be a string")
  
  if(!is.character(consumer_secret))
    stop("consumer_secret must be a string")
  
  url<-paste("http://api.tumblr.com/v2/blog/",base_hostname,"/followers",sep="")
  bodyParams<-list(limit=limit,offset=offset)
  connection<-"GET"
  
  res<-http.connection(url,token,bodyParams,consumer_key,consumer_secret,connection)
  res<-fromJSON(res)
  
  return(res)
  
}
