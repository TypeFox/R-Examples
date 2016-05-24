#' To return basic playlist
#'
#' @param api_key Echo Nest API key
#' @param type the type of the playlist to be generated
#' @param artist_id artist id 
#' @param artist artist name
#' @param song_id song ID
#' @param genre genre name
#' @param track_id track ID
#' @param results the number of results desired
#' @param partner partner catalog
#' @param tracks tracks info
#' @param limited_interactivity interactivity limitation
#' @return data frame giving basic playlist
#' @export
#' @examples
#' \dontrun{
#' data=basic_playlist(api_key,type="artist-radio",artist=c("coldplay","adele"))
#' }

basic_playlist=function(api_key,type=NA,artist_id=NA,artist=NA,song_id=NA,
                        genre=NA,track_id=NA,results=15,partner=NA,tracks=F,limited_interactivity=NA)
{
  url=paste("http://developer.echonest.com/api/v4/playlist/basic?api_key=",api_key,"&format=json",sep="")
  if(results>100)
  {
    stop("results should be less than or equal to 100")  
  }
  if(!is.na(type))
  {
    type=gsub(" ","+",type)
    url=paste(url,"&type=",type,sep="")
  }
  
  if(!is.na(artist_id))
  {
    len=length(artist_id)
    for(i in 1:len)
    {
      artist_id[i]=gsub(" ","+",artist_id[i])
      url=paste(url,"&artist_id=",artist_id[i],sep="")
      
    }
  }
  
  if(!is.na(artist))
  {
    len=length(artist)
    for(i in 1:len)
    {
      artist[i]=gsub(" ","+",artist[i])
      url=paste(url,"&artist=",artist[i],sep="")
      
    }
  }
  
  if(!is.na(song_id))
  {
    len=length(song_id)
    for(i in 1:len)
    {
      song_id[i]=gsub(" ","+",song_id[i])
      url=paste(url,"&song_id=",song_id[i],sep="")
      
    }
  }
  
  if(!is.na(genre))
  {
    len=length(genre)
    for(i in 1:len)
    {
      genre[i]=gsub(" ","+",genre[i])
      url=paste(url,"&genre=",genre[i],sep="")
      
    }
  }
  
  if(!is.na(track_id))
  {
    len=length(track_id)
    for(i in 1:len)
    {
      track_id[i]=gsub(" ","+",track_id[i])
      url=paste(url,"&track_id=",track_id[i],sep="")
      
    }
  }
  
  if(tracks)
  {
    url=paste(url,"&bucket=tracks",sep="")
  }
  
  
  url=paste(url,"&results=",results,sep="")
  
  if(!is.na(limited_interactivity))
  {
    url=paste(url,"&limited_interactivity=",limited_interactivity,sep="")
  }
  
  if(!is.na(partner))
  {
    url=paste(url,"&bucket=id:",partner,"&limit=true",sep="")  
  }
  rd=getURL(url)
  rd=fromJSON(rd)
  rd$response$songs
}