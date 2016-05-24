#' To returns a list of the overall top terms
#'
#' @param api_key Echo Nest API key
#' @param results the number of results desired
#' @return data frame giving top terms
#' @export
#' @examples
#' \dontrun{
#' data=get_top_terms(api_key)
#' }

get_top_terms=function(api_key,results=NA)
{
  url=paste("http://developer.echonest.com/api/v4/artist/top_terms?api_key=",api_key,"&format=json",sep="")
  
  if(!is.na(results))
  {
    url=paste(url,"&results=",results,sep="")
  }
  
  rd=getURL(url)
  rd=fromJSON(rd)
  
  rd$response$terms
}