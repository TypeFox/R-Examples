#' To return similar genres to a given genre
#'
#' @param api_key Echo Nest API key
#' @param genre the genre name
#' @param description genre's description
#' @param urls genre's urls
#' @param results the number of results desired
#' @param start the desired index of the first result returned
#' @return data frame giving similar genres
#' @export
#' @examples
#' \dontrun{
#' data=similar_genres(api_key,genre="rock")
#' }

similar_genres=function(api_key,genre=NA,description=T,urls=T,start=NA,results=15)
{
  url=paste("http://developer.echonest.com/api/v4/genre/similar?api_key=",api_key,"&format=json",sep="")
  if(results>100)
  {
    stop("results should be less than or equal to 100")  
  }
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
  if(!is.na(start))
  {
    url=paste(url,"&start=",start,sep="")
  }
  
  url=paste(url,"&results=",results,sep="")
  
  rd=getURL(url)
  rd=fromJSON(rd)
  rd$response$genres
}