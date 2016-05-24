#' To search song
#'
#' @param api_key Echo Nest API key
#' @param artist artist's name
#' @param artist_id artist's id
#' @param hotttnesss song's hotttnesss
#' @param combined query both artist and title fields
#' @param style artist's style
#' @param title song's title
#' @param sort to sort ascending or descending
#' @param audio_summary song's audio summary
#' @param partner partner catalog
#' @param discovery artist's discovery measure
#' @param min_name features' minimum value settings
#' @param min_val features' minimum value settings
#' @param max_name features' maximum value settings
#' @param max_val features' maximum value settings
#' @param start the desired index of the first result returned
#' @param results maximum size
#' @param artist_location artist location
#' @param mode the mode of songs
#' @param key the key of songs in the playlist
#' @param currency song currency
#' @param description song's description
#' @param rank_type For search by description, style or mood indicates whether results should be ranked by query relevance or by artist familiarity
#' @param mood a mood like happy or sad
#' @param familiarity song's familiarity
#' @param song_type controls the type of songs returned
#' @param artist_start_year_before Matches artists that have an earliest start year before the given value
#' @param artist_start_year_after Matches artists that have an earliest start year after the given value
#' @param artist_end_year_before 	Matches artists that have a latest end year before the given value
#' @param artist_end_year_after	Matches artists that have a latest end year after the given value
#' @return data frame giving artist's familiarity
#' @export
#' @examples
#' \dontrun{
#' data=search_songs(api_key,style="pop",results=31)
#' }

search_songs=function(api_key,artist=NA,artist_id=NA,title=NA,hotttnesss=T,style=NA,artist_location=T,
                      combined=NA,sort=NA,audio_summary=F,partner=NA,min_name=NA,discovery=T,
                      max_name=NA,min_val=NA,max_val=NA,start=NA,results=15,mode=NA,key=NA,currency=T,
                      description=NA,rank_type="relevance",mood=NA,familiarity=T,
                      song_type=NA,artist_start_year_before=NA,artist_start_year_after=NA,
                      artist_end_year_before=NA,artist_end_year_after=NA)
{
  url=paste("http://developer.echonest.com/api/v4/song/search?api_key=",api_key,"&format=json",sep="")
  final=""
  if(results>100)
  {
    stop("results should be less than or equal to 100")  
  }
  if(!is.na(artist))
  {
    artist=gsub(" ","+",artist)
    url=paste(url,"&artist=",artist,sep="")
  }
  
  if(!is.na(artist_id))
  {
    url=paste(url,"&artist_id=",artist_id,sep="")
  }
  
  if(!is.na(combined))
  {
    combined=gsub(" ","+",combined)
    url=paste(url,"&combined=",combined,sep="")
  }
  
  if(!is.na(description))
  {
    description=gsub(" ","+",description)
    url=paste(url,"&description=",description,sep="")
  }
  
  if(!is.na(song_type))
  {
    song_type=gsub(" ","+",song_type)
    url=paste(url,"&song_type=",song_type,sep="")
  }
  
  if(!is.na(mood))
  {
    mood=gsub(" ","+",mood)
    url=paste(url,"&mood=",mood,sep="")
  }
  
  if(!is.na(mode))
  {
    mode=gsub(" ","+",mode)
    url=paste(url,"&mode=",mode,sep="")
  }
  
  if(!is.na(key))
  {
    key=gsub(" ","+",key)
    url=paste(url,"&key=",key,sep="")
  }
  
  if(hotttnesss)
  {
    url=paste(url,"&bucket=song_hotttnesss&bucket=song_hotttnesss_rank&bucket=artist_hotttnesss&bucket=artist_hotttnesss_rank",sep="")
  }
  
  if(familiarity)
  {
    url=paste(url,"&bucket=artist_familiarity&bucket=artist_familiarity_rank",sep="")  
  }
  
  if(currency)
  {
    url=paste(url,"&bucket=song_currency&bucket=song_currency_rank",sep="")
  }
  
  if(discovery)
  {
    url=paste(url,"&bucket=artist_discovery&bucket=artist_discovery_rank",sep="")
  }
  
  if(audio_summary)
  {
    url=paste(url,"&bucket=audio_summary",sep="")
  }

  if(artist_location)
  {
    url=paste(url,"&bucket=artist_location",sep="")
  }
  
  
  if(!is.na(sort))
  {
    url=paste(url,"&sort=",sort,sep="")
  }

  if(!is.na(partner))
  {
    url=paste(url,"&bucket=id:",partner,"&limit=true",sep="")
  }

  if(!is.na(style))
  {
    url=paste(url,"&style=",style,sep="")
  }

  if(!is.na(min_name))
  {
    len=length(min_name)
    for(i in 1:len)
    {
      url=paste(url,"&",min_name[i],"=",min_val[i],sep="")
    }
  }

  if(!is.na(max_name))
  {
    len=length(max_name)
    for(i in 1:len)
    {
      
      url=paste(url,"&",max_name[i],"=",max_val[i],sep="")
    }
  }

url=paste(url,"&rank_type=",rank_type,sep="")  
  
  if(!is.na(title))
  {
    title=gsub(" ","+",title)
    url=paste(url,"&title=",title,sep="")
  }

if(!is.na(artist_start_year_before))
{
  url=paste(url,"&artist_start_year_before=",artist_start_year_before,sep="")
}

if(!is.na(artist_start_year_after))
{
  url=paste(url,"&artist_start_year_after=",artist_start_year_after,sep="")
}

if(!is.na(artist_end_year_before))
{
  url=paste(url,"&artist_end_year_before=",artist_end_year_before,sep="")
}

if(!is.na(artist_end_year_after))
{
  url=paste(url,"&artist_end_year_after=",artist_end_year_after,sep="")
}
if(!is.na(start))
{
  url=paste(url,"&start=",start,sep="")  
}

  url=paste(url,"&results=",results,sep="")  
  

rd=getURL(url)
rd=fromJSON(rd)

rd$response$songs

}
