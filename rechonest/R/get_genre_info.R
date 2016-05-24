
#' To get basic information about a genre
#'
#' @param api_key Echo Nest API key
#' @param genre the genre name
#' @param description genre's description
#' @param urls genre's urls
#' @return data frame giving basic info about a genre
#' @export
#' @examples
#' \dontrun{
#' data=get_genre_info(api_key,genre="post rock")
#' }

get_genre_info=function(api_key,genre,description=T,urls=T)
{
  url=paste("http://developer.echonest.com/api/v4/genre/profile?api_key=",api_key,"&format=json",sep="")
  genre=gsub(" ","+",genre)
  
  url=paste(url,"&name=",genre,sep="")
  
  if(description)
  {
    url=paste(url,"&bucket=description",sep="")
  }
  if(urls)
  {
    url=paste(url,"&bucket=urls",sep="")
  }
  rd=getURL(url)
  rd=fromJSON(rd)
  rd$response$genres
}
