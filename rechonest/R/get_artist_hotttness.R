
#' To get artist's hotttnesss
#'
#' @param api_key Echo Nest API key
#' @param name artist's name
#' @param id artist's id
#' @return data frame giving artist's hotttnesss
#' @export
#' @examples
#' \dontrun{
#' data=get_artist_hotttnesss(api_key,name="coldplay")
#' }

get_artist_hotttnesss=function(api_key,name=NA,id=NA)
{
  url=paste("http://developer.echonest.com/api/v4/artist/hotttnesss?api_key=",api_key,"&format=json",sep="")
  
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
  
  data=rd$response$artist
  data
}
