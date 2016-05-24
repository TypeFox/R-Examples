#' @title Observe your usage
#' @description \code{zeit_client} does not provide content per se, 
#' but lets you get information about your API usage.
#'
#' @param print if \code{TRUE} (default), the client information is printed.
#' @return a list of information about the client and API usage
#' @export
#' @examples
#' \dontrun{
#' zeit_client()
#' }
zeit_client <- function(print=TRUE) {
	# make request
  req <- zeit_get_url(path="client")
  raw <- zeit_parse(req)
  #reset <- as.POSIXct(raw$reset, origin="1970-01-01")
  
  # return
  if(print) print(raw)
  invisible(raw)
}


#' @title Search the ZEIT archive
#' @description \code{zeit_search} exposes a search for ZEIT archive items. 
#' You can set search queries, paginate, sort and partially select the fields, 
#' that should be returned. Articles, that match your query, are returned with
#' a reduced set of meta data.
#'
#' @param endpoint one of \code{author}, \code{content}, \code{department}, \code{keyword}, 
#' \code{product} or \code{series} -- provides specific search functionalities.
#' @param query the main search query; single string value or vector of strings.
#' @param fields partially select output fields, as string value or vector of strings for multiple fields.
#' @param limit limit the amount of matches to return; set to \code{10} by default.
#' @param offset offset for the list of matches; set to \code{0} by default.
#' @param sort sort search results by any of the returned \code{fields}. 
#' Vector of two (\code{c([field], [direction])}), giving field and direction keyword. 
#' Direction keywords are \code{asc} and \code{desc} for an ascending or descending 
#' sort order respectively. Multiple sort orders are accepted as \code{list} of such vectors.
#' @param print if \code{TRUE} (default) the search results are printed.
#' @return A list of matches to the query.
#' @details \emph{Endpoints}
#' 
#' The API is structured into several endpoints that provide specific functionalities:
#' \tabular{lll}{
#' \tab \code{author} \tab search all authors \cr
#' \tab \code{content} \tab search for content \cr
#' \tab \code{department} \tab search all departments \cr
#' \tab \code{keyword} \tab search all keywords \cr
#' \tab \code{product} \tab search all products \cr
#' \tab \code{series} \tab search all series
#' }
#'
#' \emph{Query syntax}
#' 
#' You can search the entire article text and all meta data simply by setting the query parameter 
#' to your search phrase. The search uses entire strings "as is". To search for multiple tokens 
#' use a vector of strings.
#'
#' All fields of an article can be queried individually by using \code{[field]:[search string]}. 
#' For example, to get articles that have the word "Kennedy" in their headline, you would 
#' search for \code{"title:Kennedy"}.
#'
#' Currently all \code{endpoint}s other than \code{content} only support simple search phrases 
#' with asterisk (\code{*}) wildcards.
#' @source \url{http://developer.zeit.de/docs/}
#' @export
#' @examples
#' \dontrun{
#' # simple content search
#' zeit_search(endpoint="content", query="bayreuth")
#' zeit_search("content", "bayreuth") # same same
#' 
#' # multiple tokens
#' zeit_search("content", c("bayreuth", "festspiele"))
#' 
#' # entire string
#' zeit_search("content", "bayreuther festspiele")
#' 
#' # field query
#' zeit_search("content", "title:bayreuth")
#'
#' # partial selection
#' zeit_search("content", "bayreuth", fields=c("title", "teaser_text"))
#'
#' # pagination
#' zeit_search("content", "bayreuth", limit=1) # just one match
#' zeit_search("content", "bayreuth", limit=1, offset=1) # just the second match
#'
#' # sorting
#' zeit_search("content", "bayreuth", 
#'   sort=c("release_date", "asc")) # sort by date
#' zeit_search("content", "bayreuth", 
#'   sort=list(c("release_date", "desc"), c("title", "asc"))) # sort by date and title
#'
#' # hide matches
#' bt.matches <- zeit_search("content", "bayreuth", print=FALSE)
#'
#' # author search
#' zeit_search(endpoint="author", query="Stefan Locke")
#' }
zeit_search <- function(endpoint, query, fields, limit=10, offset=0, sort, print=TRUE) {
	# prepare endpoint
	avail.endpoints <- c("author", "content", "department", "keyword", "product", "series")
	endpoint <- avail.endpoints[pmatch(endpoint, avail.endpoints)]
	
	# prepare query
	if(length(query)>1) query <- sapply(query, function(x) if(length(grep(" ", x, fixed=TRUE))!=0) paste0("\"", x, "\"") else paste(x))
	else if(length(grep(" ", query, fixed=TRUE))!=0) query <- paste0("\"", query, "\"")
	query <- paste0(query, collapse="+")
	
	# prepare fields
	if(endpoint == "author") avail.fields <- c("href", "id", "type", "uri", "value")
	else if(endpoint == "content") avail.fields <- c("subtitle", "uuid", "title", "href", "release_date", "uri", "snippet", "supertitle", "teaser_text", "teaser_title")
	else if(endpoint == "department") avail.fields <- c("uri", "value", "parent", "href")
	else if(endpoint == "keyword") avail.fields <- c("uri", "value", "lexical", "type", "score", "href")
	else if(endpoint == "product") avail.fields <- c("uri", "value")
	else if(endpoint == "series") avail.fields <- c("uri", "value", "href")
	if(!missing(fields)) {
		fields <- avail.fields[pmatch(fields, avail.fields)]
		if(length(fields)>1) fields <- paste0(fields, collapse=",")
	} else fields = NULL
	
	# prepare sort
	if(!missing(sort)) {
		if(is.list(sort)) {
			# check sort fields and direction
			if(!all(sapply(sort, length) == 2)) stop("Cannot resolve sort parameter. Please specify sort as vector 'c([field], [direction])' or list of such vectors.")
			sort.fields <- sapply(1:length(sort), function(x) sort[[x]][1])
			check.fields <- match(sort.fields, avail.fields)
			if(any(is.na(check.fields))) stop("Sort field(s) not available for \'", endpoint, "\': ", sort.fields[which(is.na(check.fields))])
			sort.dir <- sapply(1:length(sort), function(x) sort[[x]][2])
			check.dir <- match(sort.dir, c("asc", "desc"))
			if(any(is.na(check.dir))) stop("Sort direction(s) not available: ", sort.dir[which(is.na(check.dir))])
			sort <- paste(lapply(sort, function(x) paste(x, collapse=" ")), collapse=", ")
		} else {
			sort <- paste(sort, collapse=" ")
		}
	} else sort <- NULL
	
	# make request
  if(is.null(fields)) {
  	if(is.null(sort)) req <- zeit_get_url(path=endpoint, q=query, limit=limit, offset=offset)
  	else req <- zeit_get_url(path=endpoint, q=query, limit=limit, offset=offset, sort=sort)
  } else {
  	if(is.null(sort)) req <- zeit_get_url(path=endpoint, q=query, fields=fields, limit=limit, offset=offset)
  	else req <- zeit_get_url(path=endpoint, q=query, fields=fields, limit=limit, offset=offset, sort=sort)
  }
  raw <- zeit_parse(req)
  
  # return
  if(print) print(raw)
  invisible(raw)
}


#' @title Get detailled content from the ZEIT archive
#' @description \code{zeit_get} will get you all available metadata for a specific item.
#'
#' @param endpoint one of \code{author}, \code{content}, \code{department}, \code{keyword}, 
#' \code{product} or \code{series} -- see \code{\link{zeit_search}}.
#' @param id item id.
#' @param fields partially select output fields, as string value or vector of strings for multiple fields.
#' @param print if \code{TRUE} (default) the meta data are printed.
#' @return List of metadata items.
#' @details \emph{Endpoints}
#' 
#' The API is structured into several endpoints that provide specific functionalities:
#' \tabular{lll}{
#' \tab \code{author} \tab content by this author \cr
#' \tab \code{content} \tab get content by ID \cr
#' \tab \code{department} \tab content from this department \cr
#' \tab \code{keyword} \tab content about this keyword \cr
#' \tab \code{product} \tab content from this product \cr
#' \tab \code{series} \tab content in this series
#' }
#' @source \url{http://developer.zeit.de/docs/}
#' @export
#' @examples
#' \dontrun{
#' # get article metadata by ID
#' zeit_get("content", "3Ed7KYJOO2MXu5SQtnudQA")
#' 
#' # partial selection of output fields
#' zeit_get("content", "3Ed7KYJOO2MXu5SQtnudQA", 
#'   fields=c("title", "release_date", "href"))
#' 
#' # hide result
#' article.meta <- zeit_get("content", "3Ed7KYJOO2MXu5SQtnudQA", print=FALSE)
#' }
zeit_get <- function(endpoint, id, fields, print=TRUE) {
	# prepare endpoint
	avail.endpoints <- c("author", "content", "department", "keyword", "product", "series")
	endpoint <- avail.endpoints[pmatch(endpoint, avail.endpoints)]
	
	# prepare fields
	if(endpoint == "author") avail.fields <- c("href", "id", "type", "uri", "value")
	else if(endpoint == "content") avail.fields <- c("subtitle", "uuid", "title", "href", "release_date", "uri", "snippet", "supertitle", "teaser_text", "teaser_title")
	else if(endpoint == "department") avail.fields <- c("uri", "value", "parent", "href")
	else if(endpoint == "keyword") avail.fields <- c("uri", "value", "lexical", "type", "score", "href")
	else if(endpoint == "product") avail.fields <- c("uri", "value")
	else if(endpoint == "series") avail.fields <- c("uri", "value", "href")
	if(!missing(fields)) {
		fields <- avail.fields[pmatch(fields, avail.fields)]
		if(length(fields)>1) fields <- paste0(fields, collapse=",")
	} else fields = NULL
	
	# prepare path
	path <- paste0(endpoint, "/", id)
	
	# make request
  if(is.null(fields)) req <- zeit_get_url(path=path)
  else req <- zeit_get_url(path=path, fields=fields)
  raw <- zeit_parse(req)
  
  # return
  if(print) print(raw)
  invisible(raw)
}
