#' @title Generic extraction function which calls boilerpipe extractors
#' @description It is the actual workhorse which directly calls the boilerpipe Java library. Typically called through
#' functions as listed for parameter \code{exname}.
#' @param exname character specifying the extractor to be used. 
#' It can take one of the following values:
#' \itemize{
#' \item{\code{\link{ArticleExtractor}}}{A full-text extractor which is tuned towards news articles.}
#' \item{\code{\link{ArticleSentencesExtractor}}}{A full-text extractor which is tuned towards extracting sentences from news articles.}
#' \item{\code{\link{CanolaExtractor}}}{A full-text extractor trained on a 'krdwrd'.}
#' \item{\code{\link{DefaultExtractor}}}{A quite generic full-text extractor.}
#' \item{\code{\link{KeepEverythingExtractor}}}{Marks everything as content.}
#' \item{\code{\link{LargestContentExtractor}}}{A full-text extractor which extracts the largest text component of a page.}
#' \item{\code{\link{NumWordsRulesExtractor}}}{A quite generic full-text extractor solely based upon the number of words per block.}
#' }
#' @param content Text content or URL as character
#' @param asText should content specifed be treated as actual text to be extracted or url (from which HTML document is first downloaded and extracted afterwards), defaults to TRUE
#' @param ... additional parameters
#' @references \url{http://code.google.com/p/boilerpipe/}
#' @importFrom rJava .jnew
#' @importFrom rJava .jcall
#' @return extracted text as character
#' @author Mario Annau
#' @export 
Extractor <- function(exname, content, asText = TRUE, ...){
	
	excontent <- character(0)
	if(asText){
		excontent <- .jnew("java/lang/String", content)
	}else{ #assume that content is url
		excontent <- .jnew("java/net/URL", content)
	}
	
	expath <- paste("de/l3s/boilerpipe/extractors", exname, sep = "/")
	
	
	ex <- .jnew(expath)
	content <- .jcall(ex, returnSig = "S", method = "getText", excontent, ...)
	
	#FIXME: Encoding problems on windows workaround
	if(.Platform$OS.type == "windows"){
		content <- iconv(content, "UTF-8", "latin1")
	}
		
	content
}

#' @title A full-text extractor which is tuned towards news articles. 
#' @description In this scenario it achieves higher accuracy than \code{\link{DefaultExtractor}}.
#' @param content Text content as character
#' @param ... additional parameters
#' @seealso \code{\link{Extractor}}
#' @return extracted text as character
#' @examples
#' data(content)
#' extract <- ArticleExtractor(content)
#' @author Mario Annau
#' @export 
ArticleExtractor <- function(content, ...){
	Extractor("ArticleExtractor", content, ...)
}

#' @title A full-text extractor which is tuned towards extracting sentences from news articles.
#' @param content Text content as character
#' @param ... additional parameters
#' @seealso \code{\link{Extractor}}
#' @return extracted text as character
#' @examples
#' data(content)
#' extract <- ArticleSentencesExtractor(content)
#' @author Mario Annau
#' @export 
ArticleSentencesExtractor <- function(content, ...){
	Extractor("ArticleSentencesExtractor", content, ...)
}

#' @title A full-text extractor trained on a 'krdwrd' Canola (see \code{https://krdwrd.org/trac/attachment/wiki/Corpora/Canola/CANOLA.pdf}.
#' @param content Text content as character
#' @param ... additional parameters
#' @seealso \code{\link{Extractor}}
#' @return extracted text as character
#' @examples
#' data(content)
#' extract <- CanolaExtractor(content)
#' @author Mario Annau
#' @export 
CanolaExtractor <- function(content, ...){
	Extractor("CanolaExtractor", content, ...)
}

#' @title A quite generic full-text extractor.
#' @param content Text content as character
#' @param ... additional parameters
#' @seealso \code{\link{Extractor}}
#' @return extracted text as character
#' @examples
#' data(content)
#' extract <- DefaultExtractor(content)
#' @author Mario Annau
#' @export 
DefaultExtractor <- function(content, ...){
	Extractor("DefaultExtractor", content, ...)
}

#' @title Marks everything as content.
#' @param content Text content as character
#' @param ... additional parameters
#' @seealso \code{\link{Extractor}}
#' @return extracted text as character
#' @examples
#' data(content)
#' extract <- KeepEverythingExtractor(content)
#' @author Mario Annau
#' @export 
KeepEverythingExtractor <- function(content, ...){
	Extractor("KeepEverythingExtractor", content, ...)
}


# FIXME: Some issues with kMin Parameter
#KeepEverythingWithMinKWordsExtractor <- function(content, kMin = 20, ...){
#	kMin.integer <- .jnew("java/lang/Integer", as.integer(kMin))
#	#kMin.integer$parseInt(as.character(kMin))
#	Extractor("KeepEverythingWithMinKWordsExtractor", content, kMin = kMin.integer, ...)
#}

#' @title A full-text extractor which extracts the largest text component of a page.
#' @description For news articles, it may perform better than the \code{\link{DefaultExtractor}},
#' but usually worse than \code{\link{ArticleExtractor}}.
#' @param content Text content as character
#' @param ... additional parameters
#' @seealso \code{\link{Extractor}}
#' @return extracted text as character
#' @examples
#' data(content)
#' extract <- LargestContentExtractor(content)
#' @author Mario Annau
#' @export 
LargestContentExtractor <- function(content, ...){
	Extractor("LargestContentExtractor", content, ...)
}

#' @title A quite generic full-text extractor solely based upon the number of words per block (the current, the previous and the next block).
#' @param content Text content as character
#' @param ... additional parameters
#' @seealso \code{\link{Extractor}}
#' @return extracted text as character
#' @examples
#' data(content)
#' extract <- NumWordsRulesExtractor(content)
#' @author Mario Annau
#' @export 
NumWordsRulesExtractor <- function(content, ...){
	Extractor("NumWordsRulesExtractor", content, ...)
}