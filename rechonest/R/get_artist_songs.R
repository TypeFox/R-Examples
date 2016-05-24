
#' To get artist's songs
#'
#' @param api_key Echo Nest API key
#' @param name artist's name
#' @param id artist's id
#' @param results maximum size
#' @param start the desired index of the first result returned
#' @return data frame giving artist's songs
#' @export
#' @examples
#' \dontrun{
#' data=get_artist_songs(api_key,name="coldplay")
#' }

get_artist_songs=function(api_key,name=NA,id=NA,start=NA,results=15)
{
  url=paste("http://developer.echonest.com/api/v4/artist/songs?api_key=",api_key,"&format=json",sep="")
  final=""
  if(results>100)
  {
    stop("results should be less than or equal to 100")  
  }
  if(!is.na(name))
  {
    name=gsub(" ","+",name)
    url=paste(url,"&name=",name,sep="")
  }
  if(!is.na(id))
  {
    url=paste(url,"&id=",id,sep="")
  }
  
  if(!is.na(start))
  {
    url=paste(url,"&start=",start,sep="")
  }
  
  
  url=paste(url,"&results=",results,sep="")
  
  rd=getURL(url)
  rd=fromJSON(rd)
  
  total=rd$response$total
  data=rd$response$songs
  final=data
  
  final$total=total
  
  final=as.data.frame(final)
  final
  
  final
}
