posts <-
function(base_hostname=NA,limit=20,offset=0,api_key=NA,type=NA,id=NA,tag=NA,
                reblog_info=FALSE,notes_info=FALSE,filter="HTML"){
  
  if(!is.character(base_hostname))
    stop("base_hostname must be a string")
  
  if(!is.numeric(limit) || (limit<=0 || limit>=21) )
    stop("limit must be a numeric type beetwen 0 and 20 (inclusive")
  
  if(!is.numeric(offset) || (offset<0 || offset>=21) )
    stop("offset must be a numeric type greater or equal to limit")
  
  if(!is.character(api_key))
    stop("api_key must be a string")
  
  if(!(is.logical(reblog_info)))
    stop("reblog_info must be a boolean")
  
  if(!(is.logical(notes_info)))
    stop("reblog_info must be a boolean")
  
  filter_type<-c("HTML","text","raw")
  
  if(!(filter %in% filter_type))
    stop("Avaliable values for filter are:  HTML, text, raw")
  
  if(!is.character(filter))
    stop("Type must be a string")
  
  if(!is.na(type)){
    type_text=c("text", "quote", "link", "answer", "video", "audio", "photo", "chat")
    
    if(!(type %in% type_text))
      stop("Avaliable values for type are: text, quote, link, answer, video, audio, photo, chat")
    
    if(!is.character(type))
      stop("Type must be a string")
  }
  
  
  if(!is.na(id)){
    if(!is.numeric(id))
    stop("id must be a numeric type")
  }
  
    if(!is.na(tag)){
      if(!is.character(tag))
        stop("tag must be a string")
    }
  
  
  Params<-list(api_key=api_key,limit=limit,offset=offset,id=id,tag=tag,
                          reblog_info=reblog_info,notes_info=notes_info,filter=filter)
  
  len<-length(Params)
  s<-NULL
  for(i in 1:len){
    if(!is.na(Params[[i]][1]))
      s<-c(s,i)
  }
  
  bodyParams<-Params[s]
  
  
  if(is.na(type)){
    url<-paste("http://api.tumblr.com/v2/blog/",base_hostname,"/posts?",paste0(names(bodyParams),"=", as.character(bodyParams),collapse="&"),sep="")
  } else{
    url<-paste("http://api.tumblr.com/v2/blog/",base_hostname,"/posts/",as.character(type),"?",paste0(names(bodyParams),"=", as.character(bodyParams),collapse="&"),sep="")
  }
  
  return(content(GET(url)))
  
}
