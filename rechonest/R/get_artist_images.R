
#' To get artist's images
#'
#' @param api_key Echo Nest API key
#' @param name artist name
#' @param id Echo Nest ID
#' @param start the desired index of the first result returned
#' @param results the number of results desired
#' @param license the desired licenses of the returned images
#' @return data frame giving artist's images
#' @export
#' @examples
#' \dontrun{
#' data=list_genres(api_key)
#' }

get_artist_images=function(api_key,name=NA,id=NA,start=NA,results=15,license="unknown")
{
  url=paste("http://developer.echonest.com/api/v4/artist/images?api_key=",api_key,"&format=json",sep="")
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
  
  if(!is.na(id))
  {
    url=paste(url,"&id=",id,sep="")
  }
  
  url=paste(url,"&license=",license,sep="")
  
  if(!is.na(start))
  {
    url=paste(url,"&start=",start,sep="")  
  }
  
  url=paste(url,"&results=",results,sep="")  
  
  
  rd=getURL(url)
  rd=fromJSON(rd)
  
  rd$response$images
}  