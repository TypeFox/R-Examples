
#' To get genre's list
#'
#' @param api_key Echo Nest API key
#' @return data frame giving genre's list
#' @export
#' @examples
#' \dontrun{
#' data=list_genres(api_key)
#' }

list_genres=function(api_key)
{
  url=paste("http://developer.echonest.com/api/v4/genre/list?api_key=",api_key,"&format=json",sep="")
  rd=getURL(url)
  rd=fromJSON(rd)
  
  url=paste("http://developer.echonest.com/api/v4/genre/list?api_key=",api_key,"&format=json&results=1000&start=1000",sep="")
  rd1=getURL(url)
  rd1=fromJSON(rd1)
  
  data=rbind(rd$response$genres,rd1$response$genres)
  data
}
