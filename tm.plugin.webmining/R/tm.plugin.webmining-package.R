#' tm.plugin.webmining facilitates the retrieval of textual data through various 
#' web feed formats like XML and JSON. Also direct retrieval from HTML 
#' is supported. As most (news) feeds only incorporate small fractions
#' of the original text tm.plugin.webmining goes a step further and even
#' retrieves and extracts the text of the original text source.
#' Generally, the retrieval procedure can be described as a two--step process:
#' \describe{
#' \item{Meta Retrieval}{In a first step, all relevant meta feeds are retrieved.
#' From these feeds all relevant meta data items are extracted.
#' }
#' \item{Content Retrieval}{In a second step the relevant source content is retrieved.
#' Using the \code{boilerpipeR} package even the main content of \code{HTML} pages can
#' be extracted.
#' }}
#' 
#' @name tm.plugin.webmining-package
#' @aliases tm.plugin.webmining webmining
#' @docType package
#' @title Retrieve structured, textual data from various web sources
#' @author Mario Annau \email{mario.annau@@gmail}
#' @keywords package
#' @seealso \code{\link{WebCorpus}} \code{\link{GoogleFinanceSource}} \code{\link{GoogleNewsSource}} \code{\link{NYTimesSource}} \code{\link{ReutersNewsSource}} \code{\link{YahooFinanceSource}} \code{\link{YahooInplaySource}} \code{\link{YahooNewsSource}} 
#' @examples
#' \dontrun{
#' googlefinance <- WebCorpus(GoogleFinanceSource("NASDAQ:MSFT"))
#' googlenews <- WebCorpus(GoogleNewsSource("Microsoft"))
#' nytimes <- WebCorpus(NYTimesSource("Microsoft", appid = nytimes_appid))
#' reutersnews <- WebCorpus(ReutersNewsSource("businessNews"))
#' yahoofinance <- WebCorpus(YahooFinanceSource("MSFT"))
#' yahooinplay <- WebCorpus(YahooInplaySource())
#' yahoonews <- WebCorpus(YahooNewsSource("Microsoft"))
#' }
NULL

#' WebCorpus retrieved from Yahoo! News for the search term "Microsoft"
#' through the YahooNewsSource. Length of retrieved corpus is 20.
#' @name yahoonews
#' @docType data
#' @author Mario Annau
#' @keywords data
#' @examples
#' #Data set has been generated as follows:
#' \dontrun{
#' yahoonews <- WebCorpus(YahooNewsSource("Microsoft"))
#' }
NULL

#' AppID for the NYtimes-API.
#' 
#' USED ONLY FOR PACKAGE TESTING. PLEASE DOWNLOAD YOUR OWN KEY AT \url{http://developer.nytimes.com/}!!!
#' @name nytimes_appid
#' @docType data
#' @author Mario Annau
#' @keywords data
NULL
