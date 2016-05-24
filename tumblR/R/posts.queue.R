posts.queue <-
function(base_hostname=NA,limit=20,offset=0,filter="HTML",
                      token=NA,consumer_key=NA,consumer_secret=NA){
  
  if(!is.character(base_hostname))
    stop("base_hostname must be a string")
  
  if(!is.numeric(limit) || (limit<=0 || limit>=21) )
    stop("limit must be a numeric type beetwen 0 and 20 (inclusive")
  
  if(!is.numeric(offset) || (offset<0 || offset>=21) )
    stop("offset must be a numeric type greater or equal to limit")
  
  filter_type<-c("HTML","text","raw")
  
  if(!(filter %in% filter_type))
    stop("Avaliable values for filter are:  HTML, text, raw")
  
  if(!is.character(filter))
    stop("Type must be a string")
  
  if(class(token)[1]!="Token1.0")
    stop("token must be a Token1.0 type")
  
  if(!is.character(consumer_key))
    stop("consumer_key must be a string")
  
  if(!is.character(consumer_secret))
    stop("consumer_secret must be a string")
  
  Params<-list(limit=limit,offset=offset,filter=filter)
  
  len<-length(Params)
  s<-NULL
  for(i in 1:len){
    if(!is.na(Params[[i]][1]))
      s<-c(s,i)
  }
  
  bodyParams<-Params[s]
   
  url<-paste("http://api.tumblr.com/v2/blog/",base_hostname,"/posts/queue",sep="")
  connection<-"GET"
  
  res<-fromJSON(http.connection(url,token,bodyParams,consumer_key,consumer_secret,connection))
  
  return(res)
}
