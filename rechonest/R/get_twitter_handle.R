
#' To get the twitter handle for an artist
#'
#' @param api_key Echo Nest API key
#' @param name artist name
#' @param id Echo Nest ID
#' @return data frame giving twitter handle
#' @export
#' @examples
#' \dontrun{
#' data=get_twitter_handle(api_key,name="coldplay")
#' }

get_twitter_handle=function(api_key,name=NA,id=NA)
{
  url=paste("http://developer.echonest.com/api/v4/artist/twitter?api_key=",api_key,"&format=json",sep="")
  
  if(!is.na(name))
  {
    name=gsub(" ","+",name)
    url=paste(url,"&name=",name,sep="")
  }
  
  if(!is.na(id))
  {
    url=paste(url,"&id=",id,sep="")
  }
  
  rd=getURL(url)
  rd=fromJSON(rd)
  rd$response$artist
}
