
#' To get a list of artist biographies
#'
#' @param api_key Echo Nest API key
#' @param name artist name
#' @param id Echo Nest ID
#' @param results the number of results desired
#' @param start the desired index of the first result returned
#' @param license the desired licenses of the returned images
#' @return data frame giving artist's biographies
#' @export
#' @examples
#' \dontrun{
#' data=get_artist_biographies(api_key,name="coldplay")
#' }

get_artist_biographies=function(api_key,name=NA,id=NA,start=NA,results=15,license="unknown")
{
  url=paste("http://developer.echonest.com/api/v4/artist/biographies?api_key=",api_key,"&format=json",sep="")
  if(results>100)
  {
    stop("results should be less than or equal to 100")  
  }
  final=""
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
  
  total=rd$response$total
  data=rd$response$biographies
  final=data
  
  final$total=total
  
  final=as.data.frame(final)
  final
}