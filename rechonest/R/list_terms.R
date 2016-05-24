
#' To get a list of the best typed descriptive terms
#'
#' @param api_key Echo Nest API key
#' @param type term type
#' @return data frame giving best typed descriptive terms
#' @export
#' @examples
#' \dontrun{
#' data=list_terms(api_key)
#' }

list_terms=function(api_key,type="style")
{
  url=paste("http://developer.echonest.com/api/v4/artist/list_terms?api_key=",api_key,"&format=json",sep="")
  url=paste(url,"&type=",type,sep="")
  rd=getURL(url)
  rd=fromJSON(rd)
  rd$response$terms
}