
#' To search for genres by name
#'
#' @param api_key Echo Nest API key
#' @param genre the genre name
#' @param description genre's description
#' @param urls genre's urls
#' @param results the number of results desired
#' @return data frame giving searched genres
#' @export
#' @examples
#' \dontrun{
#' data=search_genre(api_key,genre="rock")\
#' }

search_genre=function(api_key,genre=NA,description=T,urls=T,results=15)
{
  url=paste("http://developer.echonest.com/api/v4/genre/search?api_key=",api_key,"&format=json",sep="")
  
  if(!is.na(genre))
  {
    genre=gsub(" ","+",genre)
    url=paste(url,"&name=",genre,sep="")
  }
  if(description)
  {
    url=paste(url,"&bucket=description",sep="")
  }
  if(urls)
  {
    url=paste(url,"&bucket=urls",sep="")
  }
  
  url=paste(url,"&results=",results,sep="")
  
  rd=getURL(url)
  rd=fromJSON(rd)
  rd$response$genres
}
