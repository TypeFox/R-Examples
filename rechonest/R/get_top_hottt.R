#' To return a list of the top hottt artists
#'
#' @param api_key Echo Nest API key
#' @param genre the set of genres of interest
#' @param start the desired index of the first result returned
#' @param results the number of results desired
#' @return data frame giving top hottt artists
#' @export
#' @examples
#' \dontrun{
#' data=get_top_hottt(api_key)
#' }

get_top_hottt=function(api_key,genre=NA,start=NA,results=15)
{
  url=paste("http://developer.echonest.com/api/v4/artist/top_hottt?api_key=",api_key,"&format=json",sep="")
  final=""
  if(results>100)
  {
    stop("results should be less than or equal to 100")  
  }
  if(!is.na(genre))
  {
    genre=gsub(" ","+",genre)
    url=paste(url,"&genre=",genre,sep="")
  }
  
  url=paste(url,"&bucket=hotttnesss",sep="")
  
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