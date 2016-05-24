info.blog <-
function(base_hostname=NA,api_key=NA){
  
  if(!is.character(base_hostname))
    stop("base_hostname must be a string")
    
  if(!is.character(api_key))
    stop("api_key must be a string")
 
  url<-paste("http://api.tumblr.com/v2/blog/",base_hostname,"/info?api_key=",api_key,sep="")
  
  return(fromJSON(getURL(url)))
  
}
