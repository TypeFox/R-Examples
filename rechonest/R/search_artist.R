
#' To search artist by using name
#'
#' @param api_key Echo Nest API key
#' @param name artist's name
#' @param style artist's style
#' @param hotttnesss artist's hotttnesss (Default is true)
#' @param description artist's description
#' @param results maximum size
#' @param sort to sort ascending or descending
#' @param partner partner catalog
#' @param min_familiarity minimum familiarity
#' @param max_familiarity maximum familiarity
#' @param min_hotttnesss minimum hotttnesss
#' @param max_hotttnesss maximum hotttnesss
#' @param artist_start_year_before Matches artists that have an earliest start year before the given value
#' @param artist_start_year_after Matches artists that have an earliest start year after the given value
#' @param artist_end_year_before 	Matches artists that have a latest end year before the given value
#' @param artist_end_year_after	Matches artists that have a latest end year after the given value
#' @param artist_location artist location
#' @param genre genre name
#' @param start the desired index of the first result returned
#' @param mood mood like happy or sad
#' @param rank_type For search by description, style or mood indicates whether results should be ranked by query relevance or by artist familiarity
#' @param fuzzy_match if true, a fuzzy search is performed
#' @return data frame giving artist's data
#' @export
#' @examples
#' \dontrun{
#' data=search_artist(api_key,"coldplay",sort="hotttnesss-desc",results=50)
#' }

search_artist=function(api_key,name=NA,style=NA,hotttnesss=T,
                       description=NA,start=NA,results=15,sort=NA,partner=NA,
                       artist_location=NA,genre=NA,mood=NA,
                       rank_type="relevance",fuzzy_match=F,
                       max_familiarity=NA,min_familiarity=NA,
                       max_hotttnesss=NA,min_hotttnesss=NA,
                       artist_start_year_before=NA,artist_start_year_after=NA,
                       artist_end_year_before=NA,artist_end_year_after=NA)
{
  
  url=paste("http://developer.echonest.com/api/v4/artist/search?api_key=",api_key,"&format=json",sep="")
  final=""
  if(results>100)
  {
    stop("results should be less than or equal to 100")  
  }
  #################NAME##########################
  if(!is.na(name))
  {
    name=gsub(" ","+",name)
    url=paste(url,"&name=",name,sep="")
  }
  
  #################STYLE####################
  if(!is.na(style))
  {
    style=gsub(" ","+",style)
    url=paste(url,"&style=",style,sep="")
  }
  
  if(hotttnesss)
  {
    url=paste(url,"&bucket=hotttnesss",sep="")
  }
  
  if(!is.na(artist_location))
  {
    artist_location=gsub(" ","+",artist_location)
    url=paste(url,"&artist_location=",artist_location,sep="")
  }
  
  if(!is.na(mood))
  {
    mood=gsub(" ","+",mood)
    url=paste(url,"&mood=",mood,sep="")
  }
  
  if(!is.na(genre))
  {
    genre=gsub(" ","+",genre)
    url=paste(url,"&genre=",genre,sep="")
  }
  
  url=paste(url,"&rank_type=",rank_type,sep="") 
  
  if(!is.na(max_familiarity))
  {
    url=paste(url,"&max_familiarity=",max_familiarity,sep="") 
  }
  
  if(!is.na(min_familiarity))
  {
    url=paste(url,"&min_familiarity=",min_familiarity,sep="") 
  }
  
  if(!is.na(max_hotttnesss))
  {
    url=paste(url,"&max_hotttnesss=",max_hotttnesss,sep="") 
  }
  
  if(!is.na(min_hotttnesss))
  {
    url=paste(url,"&min_hotttnesss=",min_hotttnesss,sep="") 
  }
  
  if(fuzzy_match)
  {
    url=paste(url,"&fuzzy_match=true",sep="")
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
  
  ###############DESCRIPTION###############
  if(!is.na(description))
  {
    description=gsub(" ","+",description)
    url=paste(url,"&description=",description,sep="")
  }
  
  #####################PARTNER#####################
  if(!is.na(partner))
  {
    partner=gsub(" ","+",partner)
    url=paste(url,"&bucket=id:",partner,"&limit=true",sep="")
  }
  
  #################SORT###############
  if(!is.na(sort))
  {
    url=paste(url,"&sort=",sort,sep="")
  }
  
  if(!is.na(start))
  {
    url=paste(url,"&start=",start,sep="")
  }
  
  url=paste(url,"&results=",results,sep="")
  rd=getURL(url)
  rd=fromJSON(rd)
  
  data=rd$response$artists
  final=data
  final
}
