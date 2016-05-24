posts.draft <-
function(base_hostname=NA,before_id=0,filter="HTML",
                      token=NA,consumer_key=NA,consumer_secret=NA){
  
  filter_type<-c("HTML","text","raw")
  
  if(!(filter %in% filter_type))
    stop("Avaliable values for filter are:  HTML, text, raw")
  
  if(!is.character(filter))
    stop("Type must be a string")
  
  if(!is.na(before_id)){
    if(!is.numeric(before_id))
      stop("id must be a number")
  }
  
  if(class(token)[1]!="Token1.0")
    stop("token must be a Token1.0 type")
  
  if(!is.character(consumer_key))
    stop("consumer_key must be a string")
  
  if(!is.character(consumer_secret))
    stop("consumer_secret must be a string")
  
  url<-paste("http://api.tumblr.com/v2/blog/",base_hostname,"/posts/draft",sep="")
  connection<-"GET"
  
  Params<-list(before_id=before_id,filter=filter)
  
  len<-length(Params)
  s<-NULL
  for(i in 1:len){
    if(!is.na(Params[[i]][1]))
      s<-c(s,i)
  }
  
  bodyParams<-Params[s]
  
  res<-fromJSON(http.connection(url,token,bodyParams,consumer_key,consumer_secret,connection))
  
  return(res)
}
