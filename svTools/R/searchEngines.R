### Pretty printing of search results
print.search <- function (x, detailed = TRUE, ...)
{
	## Print header
	cat("Search query for '", attr(x, "query"), "'\n", sep = "")
	cat(nrow(x), "hits with score ranging from", max(x$score), "to",
		min(x$score), ":\n\n")
	if (isTRUE(detailed)) {
		cat(paste(format(1:nrow(x), width = 3), " ", x$type, ": ", x$item, " (",
			x$page, ") - score ", x$score, sep = ""), sep = "\n")
		cat("\n")
	}
	return(invisible(x))
}

### Display a page containing one item in a search list
browse <- function (object, ...)
	NextMethod("browse")

browse.search <- function (object, item = 1, ...)
{
	if (!inherits(object, "search"))
		stop("'object' must be a 'search' object")
	url <- object[item, "url"]
	if (is.na(url) || url == "")
		stop("No URL defined for this item")
	browseURL(url)
}

### Search in packages pages on the web
searchPackage <- function (query, max = 30)
{
	## Build the string to submit to the website
	if (missing(query)) return(NULL)
	
	string <- paste("http://search.r-project.org/cgi-bin/namazu.cgi?query=",
		gsub(" ", "+", query), sep = "")
	mpp <- paste("max=", max, sep = "")
	restr <- "idxname=functions"
	lang <- "lang=firefox"
	url <- paste(string, mpp, lang, restr, sep = "&")
  
	## Download the result of the query and read it into R
	temp <- tempfile()
	on.exit(unlink(temp))
	download.file(url, temp, quiet = TRUE)
	result <- readLines(temp)
	if (any(regexpr("No document matching your query", result) > 0)) return(NULL)
  
	## Process the webpage to extract relevant information
	result <- grep("<dt>.*http.*/R/library/", result, value = TRUE)
	links <- sub('^.*<a href="([^"]+)">.+$', '\\1', result)
	result <- gsub(".*/R/library/", "", result)
	result <- sub("\\)</dt>.*$", "", result)
	result <- strsplit(result, '(/html/|\\.html">R: |</a></strong> \\(score: ?)')
	result <- as.data.frame(do.call(rbind, result), stringsAsFactors = FALSE)
### TODO: translate HTML characters into plain characters for column 3
	colnames(result) <- c("item", "page", "snippet", "score") 
	result$score <- as.numeric(result$score) 
	result <- cbind(data.frame(type = rep("package", nrow(result)),
		stringsAsFactors = FALSE), result, data.frame(url = links,
		stringsAsFactors = FALSE))
  
	## Structure the result by score
	result <- result[order(result$score, decreasing = TRUE), ]
  
	return(structure(result, class = c("search", "data.frame"), query = query,
		link = url))
}                      

### Search in the R Wiki
searchWiki <- function (query, max = 30)
{
	if (missing(query)) return(NULL)
	
	tmp <- tempfile() 
	on.exit(unlink(tmp))
  
	url <- sprintf("http://rwiki.sciviews.org/doku.php?do=search&id=%s", 
		gsub(" +", "+", query))
	download.file(url, tmp, quiet = TRUE) 

	results <- readLines(tmp)
	results <- results[regexpr("search_result", results) > 0]
	results <- tail(unlist(strsplit(results, "search_result")), -1)
  
	ids <- gsub('&amp.*', '', results) 
	ids <- gsub('^.*id=', '', ids) 
    
	hits <- as.numeric(gsub('^.*class=\\"search_cnt\\">(.*) Hits</span>.*',
		'\\1', results))
  
	snippets <- gsub('.*\\"search_snippet\\">', '', results) 
	snippets <- gsub('</?[^>]*>', '', snippets) 
	snippets <- gsub('\\[[^\\]*]\\]', '', snippets) 
	snippets <- gsub('<div.*', '', snippets)
	snippets <- gsub('hideLoadBar\\([^\\)]*\\)', '', snippets)
	n <- length(hits)
	if (n > max) {
		snippets <- snippets[1:max]
		ids <- ids[1:max]
		hits <- hits[1:max]
		n <- max
	}
	return(structure(data.frame(type = rep("wiki", n), item = rep("Rwiki", n),
		page = ids, snippet = snippets, score = hits,
		url = paste("http://rwiki.sciviews.org/doku.php?id",
		ids, sep = "="), stringsAsFactors = FALSE),
		class = c("search", "data.frame"), query = query, link = url))
}

### Search in the R Graph Gallery
searchGraph <- function (query, max = 30)
{
	if (missing(query)) return(NULL)
	
	url <- sprintf("http://addictedtor.free.fr/graphiques/simplesearch.php?q=%s", 
		gsub(" +", "+", query))
	tmp <- tempfile() 
	on.exit(unlink(tmp))
	download.file(url, tmp, quiet = TRUE) 

	results <- readLines(tmp)
	ids <- gsub(',.*', '', results) 
	titles <- gsub('^[0-9]*,', '', results)
	n = length(titles)
	if (n > max) {
		titles <- titles[1:max]
		ids <- ids[1:max]
		n <- max
	}
	return(structure(data.frame(type = rep("graph", n), item = titles,
		page = ids, snippet = rep("", n), score = rep(10, n),
		url = paste("http://addictedtor.free.fr/graphiques/RGraphGallery.php?graph",
		ids, sep = "="), stringsAsFactors = FALSE),
		class = c("search", "data.frame"), query = query, link = url))
}

### Convert R news/R journal bib database into a data.frame
.rbib <- function (url = "http://journal.r-project.org/RJournal.bib",
add.bibRNews = TRUE)
{
	if (isTRUE(add.bibRNews)) {
		bibRNews <- data.frame()
		load(system.file("data", "bibRNews.rda", package = "svTools"))
		if (is.null(url)) return(bibRNews)
	}
	
	tmp <- tempfile() 
	on.exit(unlink(tmp))
	download.file(url, tmp, quiet = TRUE) 
	biblines <- readLines(tmp)
    
	## Read the pdf urls (for Rnews.bib)
	strings <- biblines[regexpr("@String\\{", biblines) > 0]
	strings <- gsub('(@String\\{|\\}|")', '', strings)   
	x <- strsplit(strings, "[[:space:]]*=[[:space:]]*")   
	refs <- sapply(x, "[", 1)
	urls <- sapply(x, "[", 2)
	names(urls) <- refs
    
	begin <- grep("^@Article", biblines)
	end <- grep("^\\}", biblines)
	articles <- mapply(function (x, y) {
		lines <- biblines[x:y]
		eqLines <- c(grep("=", lines), length(lines) + 1)
		txt <- mapply(function (x1, x2) {  
				paste(lines[x1:x2], collapse = " ")    
			}, head(eqLines, -1), tail(eqLines - 1, -1))
		txt <- strsplit(txt, "[[:space:]]*=[[:space:]]*")
		fields <- gsub("[[:space:]]+", "", sapply(txt, "[", 1))
		content <- gsub("(^\\{|\\}?,?$)", "", sapply(txt, "[", 2))
		content <- gsub(" +", " ",content)
		names(content) <- fields
		if (any(fields == "pdf") && content["pdf"] %in% refs)
			content["url"] <- urls[content["pdf"]]
		## Keep only these fields
		content <- content[c("author", "title", "journal", "year", "volume",
			"number", "pages", "month", "url")]
		return(content)
	}, begin + 1, end - 1)
  
	out <- as.data.frame(t(articles), stringsAsFactors = FALSE) 
	out$year <- as.numeric(out$year)
	out$volume <- as.numeric(out$volume)
	out$number <- as.numeric(out$number)
	#colnames(out) <- uNames 
	out$issue <- paste("vol. ", out$volume, "/", out$number, " - ", out$month,
		" ", out$year, sep = "")
	
	## Do we add R News data too?
	if (isTRUE(add.bibRNews)) out <- rbind(out, bibRNews)
	return(out)
}

### Search in R bibliography (R News and R Journal)
searchBiblio <- function (query, max = 30, url = "http://journal.r-project.org/RJournal.bib",
	add.bibRNews = TRUE, ...)
{
	if (missing(query)) return(NULL)
	bib <- .rbib(url = url, add.bibRNews = add.bibRNews)
	matches <- agrep(query, paste(bib$title, bib$author), ignore.case = TRUE, ...)
	if (length(matches) == 0) return(NULL)
	subs <- bib[matches, ]
	n = nrow(subs)
	if (n > max) {
		subs <- subs[1:max, ]
		n <- max
	}
	return(structure(data.frame(type = rep("biblio", n),
		item = paste(subs$journal, subs$issue),
		page = subs$pages, snippet = paste(subs$title, " (by ", subs$author, ")",
		sep = ""), score = rep(20, n), url = subs$url, stringsAsFactors = FALSE),
		class = c("search", "data.frame"), query = query, link = url))
}

### Search the gmane mailing list
searchMailing <- function (query, max = 30, groups = "*", prefix = "gmane.comp.lang.r")
{  
	## Building the search url
	url <- sprintf("http://search.gmane.org/?query=%s&group=%s.%s&sort=relevance", 
		gsub(" +", "+", query), prefix, groups)  
  
	## Make the search
	tmp <- tempfile()
	on.exit(unlink(tmp))
	download.file(url, destfile = tmp, quiet = TRUE)
	rl <- readLines(tmp)
	## Data is returned 10 items at a time => retrieve the other pages too...
	morePages <- (max - 1) %/% 10
	if (morePages > 0) {
		morePages <- 2:(morePages + 1)
		for (page in morePages) {
			download.file(paste(url, page, sep = "&[="), destfile = tmp,
				quiet = TRUE)
			rl <- c(rl, readLines(tmp))
		}
	}
	firstLines <- rl[regexpr("^<A HREF", rl) > 0]
	links <- gsub('^<A HREF="([^"]+)".*', '\\1', firstLines)
	links <- gsub("/match=.*$", "", links)
	groups <- gsub(sprintf("http://article.gmane.org/%s.([^/]+)/.*", prefix), 
		"\\1", links)  
	relevances <- gsub('.*\\((.*)%\\)$', '\\1', firstLines)
	firstLines <- sub("<b[^>]+>", "", firstLines)
	firstLines <- sub("</b>", "", firstLines)
	titles <- gsub('.*>([^<]+)<.*', '\\1', firstLines)

	n = length(links)
	## If n > max, limit results
	if (n > max) {
		groups <- groups[1:max]
		titles <- titles[1:max]
		relevances <- relevances[1:max]
		links <- links[1:max]
		n <- max
	}
	return(structure(data.frame(type = rep("mailing", n),
		item = groups, page = titles, snippet = rep("", n), score = relevances,
		url = links, stringsAsFactors = FALSE),
		class = c("search", "data.frame"), query = query, link = url))
}
