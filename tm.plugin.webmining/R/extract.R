#' @title Extract main content from \code{TextDocument}s.
#' @description Use implemented extraction functions (through boilerpipeR) to extract main content from
#' \code{TextDocument}s.
#' @param x PlainTextDocument
#' @param extractor default extraction function to be used, defaults to \code{\link{extractContentDOM}}
#' @param ... additional parameters to extractor function
#' @export
#' @aliases extract.PlainTextDocument
extract <- function(x, extractor, ...) UseMethod("extract", x)


#' Extract Main Content from Text Documents
#' Use implemented extraction functions (through boilerpipeR) to extract main content from
#' \code{TextDocument}s.
#' @param x PlainTextDocument
#' @param extractor default extraction function to be used, defaults to \code{\link{extractContentDOM}}
#' @param ... additional parameters to extractor function
#' @importFrom NLP content
#' @noRd
#' @export
extract.PlainTextDocument <- function(x, extractor = extractContentDOM, ...){
	content(x) <- tryCatch(extractor(x, ...),
			error = function(e){
				warning(e)
				content(x)
			})
	x
} 

#' @title Simply strip HTML Tags from Document
#' @description \code{extractHTMLStrip} parses an url, character or filename, reads the DOM
#' tree, removes all HTML tags in the tree and outputs the source text without
#' markup.
#' @author Mario Annau
#' @param url character, url or filename
#' @param asText specifies if url parameter is a \code{character}, defaults to TRUE
#' @param encoding specifies local encoding to be used, depending on platform
#' @param ... Additional parameters for \code{\link{htmlTreeParse}} 
#' @seealso \code{\link{xmlNode}}
#' @importFrom XML htmlTreeParse toString.XMLNode xmlChildren xmlValue free
#' @seealso \code{\link{htmlTreeParse}} \code{\link{encloseHTML}}
#' @note Input text should be enclosed in <html>'TEXT'</html> tags to ensure correct
#' DOM parsing (issue especially under .Platform$os.type = 'windows')
#' @export
extractHTMLStrip <-
function(url, asText = TRUE, encoding, ...){
	if(missing(encoding)){
		encoding <- switch(.Platform$OS.type,
				unix = "UTF-8",
				windows = "latin1")
	}	

	if(url == ""){
		return("")
	}

	parseerror <- capture.output(tree <- htmlTreeParse(url, asText = asText, useInternalNodes = TRUE, encoding = encoding, ...))
	
	children <- xmlChildren(tree)
	childlen <- sapply(children, function(x) nchar(toString.XMLNode(x)))
	childidx <- max(which(childlen == max(childlen)))
	#childidx <- min(childidx, length(children))
	html <- children[[childidx]]
	val <- xmlValue(html)
	XML::free(tree)
	return(val)
}


#' @title Extract Main HTML Content from DOM
#' @description Function extracts main HTML Content using its Document Object Model.
#' Idea comes basically from the fact, that main content of an HTML Document
#' is in a subnode of the HTML DOM Tree with a high text-to-tag ratio.
#' Internally, this function also calls 
#' \code{assignValues}, \code{calcDensity}, \code{getMainText} 
#' and \code{removeTags}.
#' @author Mario Annau
#' @param url character, url or filename
#' @param threshold threshold for extraction, defaults to 0.5
#' @param asText boolean, specifies if url should be interpreted as character
#' @param ... Additional Parameters to \code{\link{htmlTreeParse}}
#' @seealso \code{\link{xmlNode}}
#' @references 	\url{http://www.elias.cn/En/ExtMainText}, 
#' 				\url{http://ai-depot.com/articles/the-easy-way-to-extract-useful-text-from-arbitrary-html/}
#' 				\cite{Gupta et al., DOM-based Content Extraction of HTML Documents},\url{http://www2003.org/cdrom/papers/refereed/p583/p583-gupta.html}
#' @importFrom XML xmlChildren
#' @importFrom XML toString.XMLNode
#' @importFrom XML htmlTreeParse
#' @aliases assignValues calcDensity removeTags getMainText
#' @export
extractContentDOM <-
function(url, threshold, asText = TRUE, ...){
		
		# FIXME: Hack because of roxygen2 bug (dot replaced by comma):
		if(missing(threshold)){
			threshold <- 0.5
		}

		if(url == ""){
			return("")
		}
		
		parseerror <- capture.output(tree <- htmlTreeParse(url, asText = asText, useInternalNodes = TRUE, ...))
		childlen <- sapply(xmlChildren(tree), function(x) nchar(toString.XMLNode(x)))
		childidx <- which(childlen == max(childlen))
		html <- xmlChildren(tree)[[childidx]]
		tags <- c("script" , "noscript", "style")
		htmlclean <- removeTags(html, tags)
		
		htmlannotated <- assignValues(htmlclean, FUN = calcDensity, threshold)
		content <- getMainText(htmlannotated, threshold)
		return(content)
}

#' Calculate density of html text to overall length of html tree text
#' @author Mario Annau
#' @param xn object of class xmlNode
#' @param annotate Specifies if \code{xn} should be annotated, defaults to TRUE
#' @seealso \code{\link{extractContentDOM}}, \code{\link{xmlNode}}
#' @importFrom XML toString.XMLNode
#' @importFrom XML xmlValue
#' @importFrom XML addAttributes
#' @noRd 
calcDensity <-
function(xn, annotate = TRUE){
	textlen <- nchar( xmlValue(xn))
	treelen <- nchar(toString.XMLNode(xn))
	dens <- textlen / treelen
	if(annotate & inherits(xn, "XMLInternalElementNode")){
		addAttributes(xn, "dens" = dens, "textlen" = textlen, "treelen" = treelen)
	}
	return(c(dens, textlen, treelen))
}

#' Assign Values as Attributes to xmlNode
#' @author Mario Annau
#' @param t object of class xmlNode
#' @param FUN Function to be executed
#' @param threshold maximum threshold needed to step down the tree, defaults to 0.5
#' @param attribname Name of used attribute, defaults to "attrib"
#' @param recursive should tree be recursively annotated?, defaults to TRUE
#' @param mintextlen minimum textlength needed to step down the tree
#' @param ... additional arguments for FUN
#' @seealso \code{\link{extractContentDOM}}, \code{\link{xmlNode}}
#' @importFrom XML xmlApply
#' @noRd 
assignValues <-
function(t, FUN, threshold, attribname = "attrib", recursive = TRUE, mintextlen = 10, ...){
	
	# FIXME: Hack because of roxygen2 bug (dot replaced by comma):
	if(missing(threshold)){
		threshold <- 0.5
	}

	dens <- xmlApply(t, FUN)
	dens <- do.call("rbind", dens)
	#dens <- as.data.frame(dens)
	
	
	if(!recursive){
		return(t)
	}
	lapply(t[(dens[,2] > mintextlen) & (dens[,1] < threshold)], assignValues, FUN, ...)
	return(t)
	
}
#' Get Main Text from Annotated HTML Tree
#' Main Text is obtained from Tree -Subnode where threshold > threshold and 
#' textlength is at maximum
#' @author Mario Annau
#' @param xml object of class xmlNode
#' @param threshold minimum threshold needed to be considered
#' @seealso \code{\link{extractContentDOM}}, \code{\link{xmlNode}}
#' @importFrom XML xpathSApply
#' @importFrom XML xmlValue
#' @noRd 
getMainText <-
function(xml, threshold){
	# FIXME: Hack because of roxygen2 bug (dot replaced by comma):
	if(missing(threshold)){
		threshold <- 0.5
	}

	textlen <- as.numeric( xpathSApply(xml, path = "//attribute::textlen"))
	dens <- as.numeric( xpathSApply(xml, path = "//attribute::dens"))
	
	textlen[dens < threshold] <- 0
	idxmaintext <- which(textlen == max(textlen))
	if(max(textlen) == 0){
		return("")
	}
	
	content <-  xpathSApply(xml, path = paste("//*[@textlen][@dens]",sep = ""))[[idxmaintext]]
	
	cleancontent <-  xmlValue(content)
	cleancontent <- trimWhiteSpaces(cleancontent)
	
	return(cleancontent)
}

#' Remove specified tags from (XML) Document Tree.
#' Tags and all of its inner content will be removed.
#' @author Mario Annau
#' @param xmldoc xmlDoc object of class xmlDoc 
#' @param tags character vector which specifies tags to remove
#' @seealso \code{\link{extractContentDOM}}
#' @export
#' @importFrom XML getNodeSet
#' @importFrom XML removeNodes
#' @noRd 
removeTags <-
function(xmldoc, tags){
	#remove scripts tags
	xquery <- paste("//", tags, sep = "", collapse = " | ")
	scripts <-  getNodeSet(xmldoc, path = xquery)
	ret <- removeNodes(scripts , free = rep(FALSE, length(scripts)))
	removeTags <- xmldoc
}


