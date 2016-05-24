#' @title Read content from WebXMLSource/WebHTMLSource/WebJSONSource. 
#' @description \code{readWeb} is a FunctionGenerator which specifies content retrieval from a \code{\link{WebSource}} 
#' content elements. Currently, it is defined for XML, HTML and JSON feeds through \code{readWebXML},
#' \code{readWebHTML} and \code{readWebJSON}. Also content parsers (\code{xml_content}, \code{json_content})
#' need to be defined.
#' @param spec specification of content reader
#' @param doc document to be parsed
#' @param parser parser function to be used
#' @param contentparser content parser function to be used, see also \code{tm:::xml_content} or \code{json_content}
#' @param freeFUN function to free memory from parsed object (actually only relevant for XML and HTML trees)
#' @return FunctionGenerator
#' @importFrom tm FunctionGenerator PlainTextDocument
#' @aliases readWebXML readWebHTML readWebJSON json_content 
#' @export
readWeb <- FunctionGenerator(function(spec, doc, parser, contentparser, freeFUN = NULL) {
			
	parser <- parser
	contentparser <- contentparser
	freeFUN <- freeFUN
	spec <- spec
	doc <- doc

	function(elem, language, id) {
		tree <- parser(elem$content)
	
		###Set Content
		content(doc) <- if ("content" %in% names(spec)){
							content <- contentparser(tree, spec[["content"]])
						}
						else{
							character(0)
						}		

		for (n in setdiff(names(spec), "content")){
				meta(doc, n) <- contentparser(tree, spec[[n]])
			}
			
			if(!is.null(freeFUN)){
				freeFUN(tree)
			}
			doc
		}
})

#' Read content from WebXMLSource
#' @param ... additional parameters to \code{\link{readWeb}}
#' @export
#' @importFrom XML xmlInternalTreeParse free
#' @noRd 
readWebXML <- function(...){
	parser <- function(x){
		#XML::xmlInternalTreeParse(x, asText = TRUE)
		parse(x, type = "XML")
	} 
	contentparser <- xml_content
	freeFUN <- free
	readWeb(parser = parser, contentparser = contentparser, freeFUN = freeFUN, ...)
}

#' Read content from WebHTMLSource
#' @param ... additional parameters to \code{\link{readWeb}}
#' @export
#' @importFrom XML htmlTreeParse free
#' @noRd 
readWebHTML <- function(...){
	#parser <- function(x) XML::htmlTreeParse(x, asText = TRUE, useInternalNodes = TRUE)
	parser <- function(x) parse(x, type = "HTML", useInternalNodes = TRUE)
	contentparser <- function(x, cspec) xml_content(x, cspec)
	freeFUN <- free
	readWeb(parser = parser, contentparser = contentparser, freeFUN = freeFUN, ...)
}

#' Read content from WebJSONSource
#' @param ... additional parameters to \code{\link{readWeb}}
#' @export
#' @noRd 
readWebJSON <- function(...){
	parser <- function(x) identity(x)
	contentparser <- function(x, cspec) json_content(x, cspec)
	freeFUN <- rm
	readWeb(parser = parser, contentparser = contentparser, freeFUN = freeFUN, ...)
}

#' Read content from XMLSource
#' @param doc list object from which content should be retrieved
#' @param spec list field name as character
#' @noRd
#' @importFrom XML xmlValue
xml_content <- function(doc, spec) {
	type <- spec[[1]]
	fun <- switch(type,
			node = XML::xmlValue,
			attribute = identity)
	
	if (identical(type, "unevaluated"))
		spec[[2]]
	else if (identical(type, "function") && is.function(spec[[2]]))
		spec[[2]](doc)
	else
		as.character(sapply(XML::getNodeSet(doc, spec[[2]]), fun))
}

#' Read content from JSONSource
#' @param doc list object from which content should be retrieved
#' @param spec list field name as character
#' @export
#' @noRd 
json_content <- 
function (doc, spec) 
{
	type <- spec[[1]]
	fun <- switch(type, field = identity, node = identity)
	if (identical(type, "unevaluated")) 
		spec[[2]]
	else if (identical(type, "function") && is.function(spec[[2]])) 
		spec[[2]](doc)
	else{
		as.character(sapply(doc[[spec[[2]]]], 
						fun))
	} 
}

#' Read content from NYTimesSource
#' @noRd
#' @export
readNYTimes <- readWebJSON(spec = list(
		author = list("field", c("byline", "original")),
		description = list("field", "snippet"),
		datetimestamp = list("function", function(node)
					strptime(node[["pub_date"]],
							format = "%Y-%m-%dT%H:%M:%SZ",
							tz = "EST")),
		heading = list("field", c("headline", "main")),
		origin = list("field", "web_url"),
		language = list("unevaluated", "en"),
		id = list("field", "_id")),
	doc = PlainTextDocument())

#' Read content from Google...Source
#' @importFrom XML getNodeSet xmlValue
#' @importFrom NLP meta<-
#' @noRd
#' @export
readGoogle <- readWebXML(spec = list(
		heading = list("node", "//title"),
		datetimestamp = list("function", function(node){
					loc <- Sys.getlocale("LC_TIME")
					Sys.setlocale("LC_TIME", "C")
					val <- sapply(getNodeSet(node, "//pubDate"), xmlValue)
					time <- strptime(val,format = "%a, %d %b %Y %H:%M:%S",tz = "GMT")
					Sys.setlocale("LC_TIME", loc)
					time
				}),
		origin = list("node", "//link"),
		description = list("function", function(node){
					val <- sapply(getNodeSet(node, "//item/description"), xmlValue)
					extractHTMLStrip(sprintf("<html>%s</html>", val), asText = T)
				}),
		id = list("node",  "//guid")),
	doc = PlainTextDocument())

#' Read content from Yahoo RSS Source
#' @importFrom XML getNodeSet xmlValue
#' @seealso \code{\link{YahooFinanceSource}}
#' @noRd
#' @export
readYahoo <- readWebXML(spec = list(
		heading = list("node", "//title"),
		datetimestamp = list("function", function(node){
					loc <- Sys.getlocale("LC_TIME")
					Sys.setlocale("LC_TIME", "C")
					val <- sapply(getNodeSet(node, "//pubDate"), xmlValue)
					time <- strptime(val,format = "%a, %d %b %Y %H:%M:%S",tz = "GMT")
					Sys.setlocale("LC_TIME", loc)
					time
				}),
		origin = list("node", "//link"),
		description = list("node", "//item/description"),
		id = list("node",  "//guid")),
	doc = PlainTextDocument())

#' Read content from Yahoo HTML Source
#' @importFrom XML getNodeSet xmlValue
#' @seealso \code{\link{YahooNewsSource}}
#' @noRd
#' @export
readYahooHTML <- readWebHTML(spec = list(
    heading = list("node", "//div[@class='compTitle']/h3[@class='title']/a"),
    datetimestamp = list("function", function(node){
                loc <- Sys.getlocale("LC_TIME")
                Sys.setlocale("LC_TIME", "C")
                val <- sapply(getNodeSet(node, "//span[@class='tri fc-2nd ml-10']"), xmlValue)
                time <- strptime(val, format = "%b %d %H:%M %p",tz = "GMT")
                Sys.setlocale("LC_TIME", loc)
                time
            }),
    origin = list("attribute", "//div[@class='compTitle']/h3[@class='title']/a/@href"),
    author = list("node", "//span[@class='cite']"),
    description = list("node", "//div[@class='compText']/p"),
    id = list("attribute", "//div[@class='compTitle']/h3[@class='title']/a/@href")),
doc = PlainTextDocument())

#' Read content from YahooInplaySource
#' @importFrom XML getNodeSet xmlValue
#' @noRd
#' @export
readYahooInplay <- readWebHTML(spec = list(
		heading = list("node", "//b[1]"),
		id = list("node", "//b[1]"),
		content = list("node", "//p"),
		datetimestamp = list("function", function(node){
					val <- unlist(getNodeSet(node, "//b[1]", fun = xmlValue))
					substr(val, 1, regexpr("\\s", val)-1)
				}),
		ticker  = list("node", "//p/b/a")),
	doc = PlainTextDocument())




#' Read content from ReutersNewsSource
#' @importFrom XML getNodeSet xmlValue
#' @noRd
#' @export
readReutersNews <- readWebXML(spec = list(
				heading = list("node", "//title"),
				datetimestamp = list("function", function(node){
							loc <- Sys.getlocale("LC_TIME")
							Sys.setlocale("LC_TIME", "C")
							val <- sapply(getNodeSet(node, "//pubDate"), xmlValue)
							time <- strptime(val,format = "%a, %d %b %Y %H:%M:%S",tz = "GMT")
							Sys.setlocale("LC_TIME", loc)
							time
						}),
				origin = list("node", "//link"),
				description = list("function", function(node){
							val <- sapply(getNodeSet(node, "//item/description"), xmlValue)
							extractHTMLStrip(sprintf("<html>%s</html>", val), asText = T)
						}),
				id = list("node",  "//guid"),
				category = list("node", "//category")),
		doc = PlainTextDocument())



