#' @title Read Web Content and respective Link Content from feedurls.
#' @description WebSource is derived from \code{\link[tm]{Source}}. In addition to calling the
#' base \code{\link[tm]{Source}} constructor function it also retrieves the specified
#' feedurls and pre--parses the content with the parser function.
#' The fields \code{$Content}, \code{$Feedurls} \code{$Parser} and \code{$CurlOpts} are finally
#' added to the \code{Source} object.
#' @author Mario Annau
#' @param feedurls urls from feeds to be retrieved
#' @param class class label to be assigned to \code{Source} object, defaults to "WebXMLSource"
#' @param reader function to be used to read content, see also \code{\link{readWeb}}
#' @param parser function to be used to split feed content into chunks, returns list of content elements
#' @param encoding specifies default encoding, defaults to 'UTF-8'
#' @param curlOpts a named list or CURLOptions object identifying the curl options for the handle. Type \code{listCurlOptions()} for all Curl options available.
#' @param postFUN function saved in WebSource object and called to retrieve full text content from feed urls 
#' @param retrieveFeedURL logical; Specify if feedurls should be downloaded first.
#' @param ... additional parameters passed to \code{WebSource} object/structure
#' @return WebSource
#' @export
#' @importFrom XML getNodeSet xmlValue
#' @importFrom RCurl curlOptions
WebSource <- function(feedurls, class = "WebXMLSource", reader, parser, encoding = "UTF-8",
		curlOpts = curlOptions(	followlocation = TRUE, 
				maxconnects = 20,
				maxredirs = 10,
				timeout = 30,
				connecttimeout = 30), postFUN = NULL, retrieveFeedURL = TRUE, ...){

	content_raw <- NULL
	if(retrieveFeedURL) {
		content_raw <- getURL(feedurls, .opts = curlOpts)
	} else {
		content_raw <- feedurls
	}
  # Filter empty content
  content_raw <- content_raw[sapply(content_raw, nchar) > 0]
  content_parsed <- unlist(lapply(content_raw, parser), recursive = FALSE)
  structure(list(encoding = encoding, length = length(content_parsed), names = NA_character_,
              position = 0, reader = reader, content = content_parsed, feedurls = feedurls,
              parser = parser, curlOpts = curlOpts, postFUN = postFUN, retrieveFeedURL = retrieveFeedURL, ...), 
            class = unique(c(class, "WebSource", "SimpleSource")))
}


#' @title Update WebXMLSource/WebHTMLSource/WebJSONSource
#' @description Typically, update is called from \code{link{corpus.update}} and refreshes \code{$Content} in 
#' Source object.
#' @param x Source object to be updated
#' @export source.update
#' @aliases source.update.WebXMLSource source.update.WebHTMLSource source.update.WebJSONSource
source.update <- function(x){
	UseMethod("source.update", x)	
}

#'update WebSource
#' @noRd
#' @export
source.update.WebXMLSource <- 
source.update.WebHTMLSource <- 
source.update.WebJSONSource <- 
function(x) {
  content_raw <- NULL
	if(x$retrieveFeedURL) {
    content_raw <- getURL(x$feedurls, .opts = x$curlOpts)
	} else {
    content_raw <- x$feedurls
  }
  # Filter empty content
  content_raw <- content_raw[sapply(content_raw, nchar) > 0]
  
	content_parsed <- unlist(lapply(content_raw, x$parser), recursive = FALSE)
	x$content <- content_parsed
	x$position <- 0
	x
}

#' @title Get feed Meta Data from Google Finance. 
#' @description Google Finance provides business and enterprise headlines for many companies. Coverage is 
#' particularly strong for US-Markets. However, only up to 20 feed items can be retrieved.
#' @author Mario Annau
#' @param query ticker symbols of companies to be searched for, see \url{http://www.google.com/finance}.
#' Please note that Google ticker symbols need to be prefixed with the exchange name, e.g. NASDAQ:MSFT
#' @param params additional query parameters
#' @param ... additional parameters to \code{\link{WebSource}}
#' @return WebXMLSource
#' @seealso \code{\link{WebSource}}
#' @export
#' @examples
#' \dontrun{
#' corpus <- Corpus(GoogleFinanceSource("NASDAQ:MSFT"))
#' }
#' @importFrom XML xmlInternalTreeParse
#' @importFrom XML xpathSApply
#' @importFrom XML getNodeSet
#' @importFrom XML xmlValue
#' @aliases readGoogle
GoogleFinanceSource <- function(query, params = 
				list( 	hl= 'en', 
						q=query, 
						ie='utf-8', 
						start = 0, 
						num = 20, 
						output='rss'),...){
	feed <- "http://www.google.com/finance/company_news"
	parser <- function(cr){
		tree <- parse(cr, type = "XML", asText = FALSE)
		xpathSApply(tree, path = "//item")
	}
	fq <- feedquery(feed, params)
  ws <- WebSource(feedurls = fq, class = "WebXMLSource", parser = parser, reader = readGoogle, 
      postFUN = getLinkContent, retrieveFeedURL = FALSE,...)
	ws
}

#' @title Get feed data from Yahoo! Finance.
#' @description Yahoo! Finance is a popular site which provides financial news and information. It is a large source
#' for historical price data as well as financial news. Using the typical Yahoo! Finance ticker 
#' news items can easily be retrieved. However, the maximum number of items is 20. 
#' @author Mario Annau
#' @param query ticker symbols of companies to be searched for, see \url{http://finance.yahoo.com/lookup}.
#' @param params, additional query parameters, see \url{http://developer.yahoo.com/rss/}
#' @param ... additional parameters to \code{\link{WebSource}}
#' @return WebXMLSource
#' @export
#' @examples
#' \dontrun{
#' corpus <- Corpus(YahooFinanceSource("MSFT"))
#' }
#' @seealso \code{\link{WebSource}}
#' @importFrom XML xmlInternalTreeParse
#' @importFrom XML xpathSApply
#' @importFrom XML getNodeSet
#' @importFrom XML xmlValue
#' @aliases readYahoo
YahooFinanceSource <- function(query, params = 
				list(	s= query, 
						region = "US",
						lang = "en-US"), ...){
	feed <- "http://feeds.finance.yahoo.com/rss/2.0/headline"
	
	fq <- feedquery(feed, params)
	parser <- function(cr){
		tree <- parse(cr, type = "XML", asText = FALSE)
		xpathSApply(tree, path = "//item")
	}
	ws <- WebSource(feedurls = fq, class = "WebXMLSource", parser = parser, reader = readYahoo, 
      postFUN = getLinkContent, retrieveFeedURL = FALSE, ...)
	ws
}

#' @title Get feed data from Google News Search \url{http://news.google.com/}
#' @description Google News Search is one of the most popular news aggregators on the web. News
#' can be retrieved for any customized user query. Up to 100 can be retrieved per 
#' request.
#' @author Mario Annau
#' @param query Google News Search query
#' @param params, additional query parameters
#' @param ... additional parameters to \code{\link{WebSource}}
#' @return WebXMLSource
#' @seealso \code{\link{WebSource}}
#' @export
#' @examples
#' \dontrun{
#' corpus <- Corpus(GoogleNewsSource("Microsoft"))
#' }
#' @importFrom XML xmlInternalTreeParse xpathSApply getNodeSet xmlValue newXMLNamespace
GoogleNewsSource <- function(query, params = 
				list(	hl= 'en', 
						q = query, 
						ie='utf-8', 
						num = 100, 
						output='rss'), ...){
	feed <- "http://news.google.com/news"
	fq <- feedquery(feed, params)
	parser <- function(cr){
		tree <- parse(cr, type = "XML", asText = FALSE)
		nodes <- xpathSApply(tree, path = "//item")
		xmlns1 <- lapply(nodes, newXMLNamespace, "http://purl.org/dc/elements/1.1/", "dc")
		nodes
	}
	ws <- WebSource(feedurls = fq, class = "WebXMLSource", parser = parser, reader = readGoogle,
      postFUN = getLinkContent, retrieveFeedURL = FALSE, ...)
	ws
}

#' @title Get feed data from Reuters News RSS feed channels. Reuters provides numerous feed 
#' @description channels (\url{http://www.reuters.com/tools/rss}) which can be retrieved through RSS 
#' feeds. Only up to 25 items can be retrieved---therefore an alternative retrieval
#' through the Google Reader API (\code{link{GoogleReaderSource}}) could be considered.
#' @author Mario Annau
#' @param query Reuters News RSS Feed, see \url{http://www.reuters.com/tools/rss} for a list of all feeds provided. Note that only string after 'http://feeds.reuters.com/reuters/' must be given. Defaults to 'businessNews'.
#' @param ... additional parameters to \code{\link{WebSource}}
#' @return WebXMLSource
#' @seealso \code{\link{WebSource}}
#' @export
#' @examples
#' \dontrun{
#' corpus <- Corpus(ReutersNewsSource("businessNews"))
#' }
#' @importFrom XML xmlInternalTreeParse xpathSApply getNodeSet xmlValue newXMLNamespace
#' @aliases readReutersNews
ReutersNewsSource <- function(query = 'businessNews', ...){
	feed <- "http://feeds.reuters.com/reuters"
	
	fq <- paste(feed, query, sep = "/")
	parser <- function(cr){
		tree <- parse(cr, type = "XML")
		nodes <- xpathSApply(tree, path = "//item")
		xmlns1 <- lapply(nodes, newXMLNamespace, "http://rssnamespace.org/feedburner/ext/1.0", "feedburner")
		nodes
	}

	ws <- WebSource(feedurls = fq, class = "WebXMLSource", parser = parser, reader = readReutersNews, 
      postFUN = getLinkContent, ...)
	ws
}

#' @title Get news data from Yahoo! News (\url{https://news.search.yahoo.com/search/}).
#' @description Currently, only a maximum of 10 items can be retrieved.
#' @author Mario Annau
#' @param query words to be searched in Yahoo News, multiple words must be separated by '+'
#' @param params, additional query parameters, see \url{http://developer.yahoo.com/rss/}
#' @param ... additional parameters to \code{\link{WebSource}}
#' @return WebXMLSource
#' @export
#' @examples
#' \dontrun{
#' corpus <- Corpus(YahooNewsSource("Microsoft"))
#' }
#' @seealso \code{\link{WebSource}}
#' @importFrom XML xmlInternalTreeParse
#' @importFrom XML xpathSApply
#' @importFrom XML getNodeSet
#' @importFrom XML xmlValue
#' @aliases readYahooHTML
YahooNewsSource <- function(query, params = 
				list(	p= query), ...){
	feed <- "https://news.search.yahoo.com/search"
	fq <- feedquery(feed, params)
	parser <- function(cr){
		tree <- parse(cr, type = "HTML", useInternalNodes = TRUE)
		xpathSApply(tree, path = "//div[contains(@class, 'dd algo')]")
	}
	ws <- WebSource(feedurls = fq, class = "WebXMLSource", parser = parser, reader = readYahooHTML, 
      postFUN = getLinkContent, ...)
	ws
}


#' @title Get feed data from NYTimes Article Search (\url{http://developer.nytimes.com/docs/read/article_search_api_v2}). 
#' @description Excerpt from the website: "With the NYTimes Article Search API, you can search New York Times articles 
#' from 1981 to today, retrieving headlines, abstracts, lead paragraphs, links to associated multimedia 
#' and other article metadata. Along with standard keyword searching, the API also offers faceted searching. 
#' The available facets include Times-specific fields such as sections, taxonomic classifiers and controlled 
#' vocabulary terms (names of people, organizations and geographic locations)."
#' Feed retrieval is limited to 1000 items (or 100 pages).
#' @author Mario Annau
#' @param query character specifying query to be used to search NYTimes articles
#' @param n number of items, defaults to 100
#' @param count number of results per page, defaults to 10
#' @param sleep integer; Seconds to sleep between feed retrieval.
#' @param curlOpts CURLOptions; RCurl options used for feed retrieval.
#' @param appid Developer App id to be used, obtained from \url{http://developer.nytimes.com/}
#' @param params additional query parameters, specified as list, see \url{http://developer.nytimes.com/docs/read/article_search_api}
#' @param ... additional parameters to \code{\link{WebSource}}
#' @seealso \code{\link{WebSource}}, \code{\link{readNYTimes}} 
#' @export
#' @examples
#' \dontrun{
#' #nytimes_appid needs to be specified
#' corpus <- WebCorpus(NYTimesSource("Microsoft", appid = nytimes_appid))
#' }
#' @export
#' @importFrom RJSONIO fromJSON
#' @importFrom boilerpipeR ArticleExtractor
#' @aliases readNYTimes
NYTimesSource <- function(query, n = 100, appid, count = 10, 
        sleep = 1, params = 
		list(	format="json",
				q = query,
				page = 1:ceiling(n/count),
				"api-key" = appid), 
    curlOpts = curlOptions(	followlocation = TRUE, 
        maxconnects = 10,
        maxredirs = 10,
        timeout = 30,
        connecttimeout = 30), ...){
	feed <- "http://api.nytimes.com/svc/search/v2/articlesearch.json"
	fq <- feedquery(feed, params)
	
	parser <- function(cr){
		json <- parse(cr, type = "JSON")
		json$response$docs
	}
  
  start <- seq(1, length(fq), by = count)
  end <- seq(count, length(fq), by = count)
  
  feedcontent <- sapply(1:length(start), function(i) {
              fcontent <- getURL(fq[start[i]:end[i]], .opts = curlOpts)
              Sys.sleep(sleep)
              fcontent
          })
  
	ws <- WebSource(feedurls = feedcontent, class = "WebJSONSource", parser = parser, reader = readNYTimes, 
      postFUN = getLinkContent, retrieveFeedURL = FALSE, ...)
	
	ws
}

#' @title Get News from Yahoo Inplay.
#' @description Yahoo Inplay lists a range of company news provided by Briefing.com. Since Yahoo Inplay
#' does not provide a structured XML news feed, content is parsed directly from the HTML page.
#' Therefore, no further Source parameters can be specified. The number of feed items per 
#' request can vary substantially.  
#' @author Mario Annau
#' @param ... additional parameters to \code{\link{WebSource}}
#' @return WebHTMLSource
#' @export
#' @examples
#' \dontrun{
#' corpus <- Corpus(YahooInplaySource())
#' }
#' @importFrom XML htmlTreeParse
#' @importFrom XML xpathSApply
#' @aliases readYahooInplay
YahooInplaySource <- function(...){
	url <- "http://finance.yahoo.com/marketupdate/inplay"
	parser <- function(cr){
		tree <- parse(cr, useInternalNodes = T, type = "HTML")
		xp_expr = "//div[@class= 'body yom-art-content clearfix']/p"
		paragraphs = xpathSApply(tree, xp_expr)
	}
	
	ws <- WebSource(feedurls = url, class = "WebHTMLSource", parser = parser, reader = readYahooInplay, ...)
	ws
}

#' @importFrom XML saveXML
#' @noRd
#' @export
getElem.WebXMLSource <- 
getElem.WebHTMLSource <- function(x) {
	list(content = saveXML(x$content[[x$position]]), linkcontent = NULL, uri = NULL)
}

#' @noRd
#' @export
getElem.WebJSONSource <- function(x) {
	list(content = x$content[[x$position]], linkcontent = NULL, uri = NULL)
}
