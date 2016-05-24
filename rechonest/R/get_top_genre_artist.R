

#' To Return the top artists for the given genre
#'
#' @param api_key Echo Nest API key
#' @param genre the genre name
#' @return data frame top artist of the given genre
#' @export
#' @examples
#' \dontrun{
#' data=get_top_genre_artists(api_key,genre="pop")
#' }

get_top_genre_artists=function(api_key,genre)
{
  url=paste("http://developer.echonest.com/api/v4/genre/artists?api_key=",api_key,"&format=json",sep="")
  genre=gsub(" ","+",genre)
  
  url=paste(url,"&name=",genre,sep="")
  url=paste(url,"&bucket=hotttnesss",sep="")
  
  rd=getURL(url)
  rd=fromJSON(rd)
  rd$response$artists 
}
