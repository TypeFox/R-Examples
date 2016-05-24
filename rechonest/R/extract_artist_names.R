#' To extract artist names from text.
#'
#' @param api_key Echo Nest API key
#' @param text text that contains artist names
#' @param min_hotttnesss the minimum hotttnesss for returned artists
#' @param max_hotttnesss the maximum hotttnesss for returned artists
#' @param min_familiarity the minimum familiarity for returned artists
#' @param max_familiarity the maximum familiarity for returned artists
#' @param sort specified the sort order of the results
#' @param results the number of results desired
#' @return data frame giving artist's names
#' @export
#' @examples
#' \dontrun{
#' data=extract_artist_names(api_key,text="I like adele and Maroon 5")
#' }

extract_artist_names=function(api_key,text,min_hotttnesss=NA,
                              max_hotttnesss=NA,min_familiarity=NA,
                              max_familiarity=NA,sort=NA,results=NA)
{
  text=gsub(" ","+",text)
  url=paste("http://developer.echonest.com/api/v4/artist/extract?api_key=",api_key,"&format=json&text=",text,sep="")
  final=""
  
  if(!is.na(sort))
  {
    url=paste(url,"&sort=",sort,sep="")
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
  
  if(is.na(results))
  {
    rd=getURL(url)
    rd=fromJSON(rd)
    
    data=rd$response$artists
    final=data
  }
  if(!is.na(results))
  {
    url=paste(url,"&results=",results,sep="")
    rd=getURL(url)
    rd=fromJSON(rd)
    
    data=rd$response$artists
    final=data 
  }
  
  final
}