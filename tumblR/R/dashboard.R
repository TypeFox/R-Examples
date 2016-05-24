dashboard <-
function(limit=20,offset=0,type=NA,since_id=0,reblog_info=FALSE,notes_info=FALSE,
                    token=NA,consumer_key=NA,consumer_secret=NA){
  
  if(!is.numeric(limit) || (limit<=0 || limit>=21) )
    stop("limit must be a numeric type beetwen 0 and 20 (inclusive)")
  
  if(!is.numeric(offset) || (offset<0 || offset>=21) )
    stop("offset must be a numeric type greater or equal to limit")
  
  if(!is.numeric(since_id))
    stop("since_id must be a number")
  
  if(!(is.logical(reblog_info)))
    stop("reblog_info must be a boolean")
  
  if(!(is.logical(notes_info)))
    stop("reblog_info must be a boolean")
  
  if(class(token)[1]!="Token1.0")
    stop("token must be a Token1.0 type")
  
  if(!is.character(consumer_key))
    stop("consumer_key must be a string")
  
  if(!is.character(consumer_secret))
    stop("consumer_secret must be a string")
  
   
  if(!is.na(type)){
    
    format_type<-c("HTML","text","raw")
    
    if(!(type %in% format_type))
      stop("Avaliable values for filter are:  HTML, text, raw")
    
  }
  
  url<-"http://api.tumblr.com/v2/user/dashboard"
  connection<-"GET"
  
  Params<-list(limit=limit,offset=offset,type=type,since_id=since_id,
                     reblog_info=reblog_info,notes_info=notes_info)
  
  len<-length(Params)
  s<-NULL
  for(i in 1:len){
    if(!is.na(Params[[i]][1]))
      s<-c(s,i)
  }
  
  bodyParams<-Params[s]
  
  res<-toJSON(http.connection(url,token,bodyParams,consumer_key,consumer_secret,connection))
  
  return(res)
}
