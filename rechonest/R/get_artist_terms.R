
#' To get artist's terms
#'
#' @param api_key Echo Nest API key
#' @param name artist's name
#' @param id artist's id
#' @return data frame giving artist's terms
#' @export
#' @examples
#' \dontrun{
#' data=get_artist_terms(api_key,name="coldplay")
#' }

get_artist_terms=function(api_key,name=NA,id=NA)
{
  url=paste("http://developer.echonest.com/api/v4/artist/terms?api_key=",api_key,"&format=json",sep="")
  
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
  
  rd$response$terms
  
}
