tagged <-
function(api_key=NA,tag=NA,before=as.integer(Sys.time()),limit=20,filter="HTML"){
  
  if(!is.character(api_key))
    stop("api_key must be a string")
  
  if(!is.numeric(limit) || (limit<=0 || limit>=21) )
    stop("limit must be a numeric type beetwen 0 and 20 (inclusive")
  
  if(!is.integer(before))
    stop("before must be an integer")
  
  filter_type<-c("HTML","text","raw")
  
  if(!(filter %in% filter_type))
    stop("Avaliable values for filter are:  HTML, text, raw")
  
  if(is.na(tag)){
    stop("tag is a requested field")
    
  } else {
    
    if(!is.character(tag))
      stop("tag must be a string")
    
    url<-paste("http://api.tumblr.com/v2/tagged?api_key=",api_key,"&tag=",tag,"&before=",
               as.character(before),"&limit=",as.character(limit),"&filter",filter,sep="")
    
  }
    
  res<-content(GET(url))
  return(res)
}
