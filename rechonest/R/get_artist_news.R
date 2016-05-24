
#' To get news about artist
#'
#' @param api_key Echo Nest API key
#' @param name artist's name
#' @param id artist's id
#' @param results maximum size
#' @param start the desired index of the first result returned 
#' @param high_relevance 	if true only items that are highly relevant for this artist will be returned
#' @return data frame giving news about artist
#' @export
#' @examples
#' \dontrun{
#' data=get_artist_news(api_key,name="coldplay",results=35)
#' }

get_artist_news=function(api_key,name=NA,id=NA,start=NA,results=15,high_relevance=F)
{
  url=paste("http://developer.echonest.com/api/v4/artist/news?api_key=",api_key,"&format=json",sep="")
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
  
  if(!high_relevance)
  {
    url=paste(url,"&high_relevance=false",sep="")
  }
  else
  {
    url=paste(url,"&high_relevance=true",sep="")
  }
  
  if(!is.na(start))
  {
    url=paste(url,"&start=",start,sep="")
  }
  
  
  url=paste(url,"&results=",results,sep="")
  
  rd=getURL(url)
  rd=fromJSON(rd)
  
  total=rd$response$total
  data=rd$response$news
  final=data
  final$total=total
  final=as.data.frame(final)
  final
}
