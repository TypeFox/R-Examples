
#' To return standard static playlist
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
#' @param song_selection to determine how songs are selected from each artist in artist-type playlists
#' @param variety the maximum variety of artists to be represented in the playlist
#' @param distribution controls the distribution of artists in the playlist
#' @param adventurousness controls the trade between known music and unknown music
#' @param seed_catalog ID of seed catalog for the playlist
#' @param sort sorting parameter
#' @param song_type controls the type of songs returned
#' @return data frame giving standard static playlist
#' @export
#' @examples
#' \dontrun{
#' data= standard_static_playlist(api_key,type="artist-radio",artist=c("coldplay","adele"))
#' }


standard_static_playlist=function(api_key,type=NA,artist_id=NA,artist=NA,song_id=NA,
                                  genre=NA,track_id=NA,results=15,partner=NA,tracks=F,
                                  limited_interactivity=NA,song_selection=NA,variety=NA,
                                  distribution=NA,adventurousness=NA,seed_catalog=NA,sort=NA,
                                  song_type=NA)
{
  url=paste("http://developer.echonest.com/api/v4/playlist/static?api_key=",api_key,"&format=json",sep="")
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
  
  if(!is.na(song_type))
  {
    url=paste(url,"&song_type=",song_type,sep="")
  }
  
  if(!is.na(song_selection))
  {
    url=paste(url,"&song_selection=",song_selection,sep="")
  }
  
  if(!is.na(variety))
  {
    url=paste(url,"&variety=",variety,sep="")
  }
  
  if(!is.na(adventurousness))
  {
    url=paste(url,"&adventurousness=",adventurousness,sep="")
  }
  
  if(!is.na(seed_catalog))
  {
    url=paste(url,"&seed_catalog=",seed_catalog,sep="")
  }
  
  if(!is.na(distribution))
  {
    url=paste(url,"&distribution=",distribution,sep="")
  }
  
  
  rd=getURL(url)
  rd=fromJSON(rd)
  rd$response$songs
}