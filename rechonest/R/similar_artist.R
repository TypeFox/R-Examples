
#' To search similar artists by using names or IDs
#'
#' @param api_key Echo Nest API key
#' @param name artists' name (maximum upto 5 names)
#' @param id Echo Nest IDs (maximum upto 5 IDs)
#' @param hotttnesss artist's hotttnesss
#' @param results maximum size
#' @param start the desired index of the first result returned
#' @param min_familiarity minimum familiarity
#' @param max_familiarity maximum familiarity
#' @param min_hotttnesss minimum hotttnesss
#' @param max_hotttnesss maximum hotttnesss
#' @param seed_catalog seed catalog
#' @param artist_start_year_before Matches artists that have an earliest start year before the given value
#' @param artist_start_year_after Matches artists that have an earliest start year after the given value
#' @param artist_end_year_before 	Matches artists that have a latest end year before the given value
#' @param artist_end_year_after	Matches artists that have a latest end year after the given value
#' @return data frame giving similar artists' data
#' @export
#' @examples
#' \dontrun{
#' data=similar_artists(api_key,name=c("coldplay","adele","maroon 5"),results=35 )
#' }

similar_artists=function(api_key,name=NA,id=NA,seed_catalog=NA,hotttnesss=T,start=0,
                         results=15, max_familiarity=NA,min_familiarity=NA,
                         max_hotttnesss=NA,min_hotttnesss=NA,
                         artist_start_year_before=NA,artist_start_year_after=NA,
                         artist_end_year_before=NA,artist_end_year_after=NA)
{
  url=paste("http://developer.echonest.com/api/v4/artist/similar?api_key=",api_key,"&format=json",sep="")
  final=""
  
  if(results>100)
  {
    stop("results should be less than or equal to 100")  
  }
  if(results+start >100)
  {
    stop("results + start should be less than or equal to 100")  
  }
  
  
  if(!is.na(name))
  {
    len=length(name)
    for(i in 1:len)
    {
      name[i]=gsub(" ","+",name[i])
      url=paste(url,"&name=",name[i],sep="")
      
    }
  }
  
  if(!is.na(id))
  {
    len=length(id)
    for(i in 1:len)
    {
      id[i]=gsub(" ","+",id[i])
      url=paste(url,"&id=",id[i],sep="")
      
    }
  }
  
  if(!is.na(seed_catalog))
  {
    len=length(seed_catalog)
    for(i in 1:len)
    {
      seed_catalog[i]=gsub(" ","+",seed_catalog[i])
      url=paste(url,"&id=",seed_catalog[i],sep="")
      
    }
    url=paste(url,"&limit=true",sep="")
  }
  
  if(hotttnesss)
  {
    url=paste(url,"&bucket=hotttnesss",sep="")
  }
  
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
  
  data=rd$response$artists
  final=data
  final
}
