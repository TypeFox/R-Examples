likes <-
function(base_hostname=NA,limit=20,offset=0,api_key=NA){
  if(!is.character(base_hostname))
    stop("base_hostname must be a string")
  
  if(!is.numeric(limit) || (limit<=0 || limit>=21) )
    stop("limit must be a numeric type beetwen 0 and 20 (inclusive")
  
  if(!is.numeric(offset) || (offset<0 || offset>=21) )
    stop("offset must be a numeric type greater or equal to limit")
  
  if(!is.character(api_key))
    stop("api_key must be a string")
  
  url<-paste("http://api.tumblr.com/v2/blog/",base_hostname,"/likes?api_key=",api_key,
             "&limit=",as.character(limit),"&offset=",as.character(offset),sep="")
  
  return(content(GET(url)))
}
