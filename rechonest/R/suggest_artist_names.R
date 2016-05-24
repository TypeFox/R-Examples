
#' To suggest artists based upon partial names
#'
#' @param api_key Echo Nest API key
#' @param name a partial artist name
#' @param results the number of results desired (maximum 15)
#' @return data frame giving artist's names 
#' @export
#' @examples
#' \dontrun{
#' data=suggest_artist_names(api_key,"cold")
#' }

suggest_artist_names =function(api_key,name,results=NA)
{
  url=paste("http://developer.echonest.com/api/v4/artist/suggest?api_key=",api_key,"&format=json",sep="")
  name=gsub(" ","+",name)
  url=paste(url,"&name=",name,sep="")
  
  if(!is.na(results))
  {
    url=paste(url,"&results=",results,sep="")
  }
  rd=getURL(url)
  rd=fromJSON(rd)
  rd$response$artists
  
}
