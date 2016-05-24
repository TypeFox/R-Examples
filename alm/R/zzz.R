#' Concatenate author names, if present, used in other functions.
#'
#' @param x a single list element with PLoS API returned nested elements
#' @return data.frame of results, with authors concatenated to single vector.
#' @export
#' @keywords internal
concat_todf <- function(x){
	if(inherits(x$author_display, "character")){
		if(length(x$author_display) > 1){x$author_display<-paste(x$author_display, collapse=", ")} else
		{x$author_display<-x$author_display}
		data.frame(x)
	} else
		if(inherits(x$author, "character")){
			if(length(x$author) > 1){x$author<-paste(x$author, collapse=", ")} else
			{x$author<-x$author}
			data.frame(x)
		} else
	data.frame(x)
}

#' Adds elements in a list that are missing because they were not returned
#' in the PLoS API call.
#'
#' @importFrom plyr laply
#' @param x A list with PLoS API returned nested elements
#' @return A list with the missing element added with an
#' 		"na", if it is missing.
#' @export
#' @keywords internal
addmissing <- function(x){
	names_ <- names(x[[which.max(laply(x, length))]])

	bbb <- function(x){
		if(identical(names_[!names_ %in% names(x)], character(0))){x} else
		{
			xx <- rep("na", length(names_[!names_ %in% names(x)]))
			names(xx) <- names_[!names_ %in% names(x)]
			c(x, xx)
		}
	}
	lapply(x, bbb)
}

#' Function to get an API key.
#'
#' Checks first to get key from your .Rprofile file for an API key with the name
#' 		'PlosApiKey'. If it is not found, the default guest key is used.
#'
#' @param x An API key, defaults to NULL.
#' @examples \dontrun{
#' getkey()
#' }
#' @keywords internal
#' @export
getkey <- function(x = NULL) {
	if(is.null(x)){
	  key <- Sys.getenv("ALM_KEY", "")
	  if(key == "") key <- getOption("PlosApiKey")

		if(is.null(key)){
			key <- "rkfDr76z75benY3pytM1"
			message("Using default key: Please get your own API key at http://api.plos.org/")
		} else
			if(class(key)=="character"){key <- key} else
				{ stop("check your key input - it should be a character string") }
	} else {
    key <- x
  }
	key
}

#' Capitalize the first letter of a character string.
#'
#' @param s A character string
#' @param strict Should the algorithm be strict about capitalizing. Defaults to FALSE.
#' @param onlyfirst Capitalize only first word, lowercase all others. Useful for
#' 		taxonomic names.
#' @examples
#' alm_capwords(c("using AIC for model selection"))
#' alm_capwords(c("using AIC for model selection"), strict=TRUE)
#' @export
#' @keywords internal
alm_capwords <- function(s, strict = FALSE, onlyfirst = FALSE) {
	cap <- function(s) paste(toupper(substring(s,1,1)),
		{s <- substring(s,2); if(strict) tolower(s) else s}, sep = "", collapse = " " )
	if(!onlyfirst){
		sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
	} else
		{
			sapply(s, function(x)
				paste(toupper(substring(x,1,1)),
							tolower(substring(x,2)),
							sep="", collapse=" "), USE.NAMES=F)
		}
}

#' Remove null elements in a list
#' @param x A list
#' @export
#' @keywords internal
almcompact <- function (x) Filter(Negate(is.null), x)

#' List the possible alert classes
#' @export
alert_classes <- function() alert_classes_strings

alert_classes_strings <- c('Net::HTTPUnauthorized','Net::HTTPRequestTimeOut','Delayed::WorkerTimeout','DelayedJobError','Net::HTTPConflict','Net::HTTPServiceUnavailable','Faraday::ResourceNotFound','ActiveRecord::RecordInvalid','TooManyErrorsBySourceError','SourceInactiveError','TooManyWorkersError','EventCountDecreasingError','EventCountIncreasingTooFastError','ApiResponseTooSlowError','HtmlRatioTooHighError','ArticleNotUpdatedError','SourceNotUpdatedError','CitationMilestoneAlert')

alm_GET <- function(x, y, sleep=0, ...){
  Sys.sleep(time = sleep)
  out <- GET(x, query=y, ...)
  stop_for_status(out)
  tt <- content(out, as = "text")
  jsonlite::fromJSON(tt, FALSE)
}

alm_POST <- function(x, y, sleep=0, ...){
  Sys.sleep(time = sleep)
  out <- POST(x, body=y, add_headers("X-HTTP-Method-Override" = "GET"), ...)
  stop_for_status(out)
  tt <- content(out, as = "text")
  jsonlite::fromJSON(tt, FALSE)
}
