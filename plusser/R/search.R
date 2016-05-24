##' Searching in Google+ Profiles
##' 
##' This function uses the Google+ API to search for a text string in profiles. 
##' Optionally, profiles can be restricted to a certain language.
##' 
##' The number of rows of the data frame returned is somewhat ambiguous.
##' Specifying the \code{results} argument will try to get that many results.
##' But there may be less (because Google could not find more) or more (because
##' Google is organizing results on pages and it would be a waste to discard
##' them automatically). If you really depend on getting not more rows than you
##' expected, use standard selection (i.e. \code{[}) to trim the results.
##' 
##' @param q The query string to search. The string is URL encoded
##'   automatically.
##' @param language A language code. See 
##'   \url{https://developers.google.com/+/api/search#available-languages}.
##' @param results The approximate number of results that will be retrieved from
##'   Google+.
##' @param nextToken,cr These are used internally to retrieve additional pages
##'   of answers from the Google+ API. Users won't need to set these arguments.
##' @return A data frame with the user ID, display names and profile type of the profiles that
##'   met the search criteria.
##' @export
##' @seealso Google+ API documentation:
##'   \url{https://developers.google.com/+/api/latest/people/search}.
##' @examples
##' \dontrun{
##' searchProfile("cats")
##' }
searchProfile <- function(q, language=NULL, results=1, nextToken=NULL, cr=1) {
  if (is.null(get("apikey", envir=gp))) stop("Set the Google+ API key first using setAPIkey().")
  if (is.null(language)) {
    languageString <- NULL
  } else {
    languageString <- paste0("&language=", language)
  }
  this.url <- paste0(base.url,
                     "people?query=",
                     curlEscape(q),
                     languageString,
                     "&maxResults=50",
                     nextToken,
                     "&key=",
                     get("apikey", envir=gp))
  
  this.res <- fromJSON(getURL(this.url), asText=TRUE)
  this.ppl <- t(sapply(this.res[["items"]], function(x) data.frame(id=x$id, dn=x$displayName, ty=x$objectType, stringsAsFactors=FALSE)))
  cr <- cr + nrow(this.ppl)
  if (!is.null(this.res[["nextPageToken"]]) & cr < results) {
    this.nextToken <- paste0("&pageToken=", this.res[["nextPageToken"]])
    this.ppl <- rbind(this.ppl,
                  searchProfile(q, language, results, this.nextToken, cr))
  }
  return(this.ppl)
}


##' Searching for Google+ Posts
##' 
##' This function uses the Google+ API to search for a text string in posts. 
##' Optionally, search results can be restricted to a certain language.
##' 
##' The result is either a simple list of items from the page that can be parsed
##' using \code{\link{parsePost}} or a data frame with that function already
##' applied.
##' 
##' The length of the list or the number of rows of the data frame are somewhat 
##' ambiguous. Specifying the \code{results} argument will try to get that many
##' results. But there may be less (because Google could not find more) or more 
##' (because Google is organizing results on pages and it would be a waste to 
##' discard them automatically). If you really depend on getting not more rows 
##' than you expected, use standard selection (i.e. \code{[}) to trim the
##' results.
##' 
##' @param q The query string to search. The string is URL encoded
##'   automatically.
##' @param ret A string specifying the kind of return value. Either a 
##'   \code{list} of the retrieved items on the page, or that list parsed into a
##'   \code{data.frame}.
##' @param language A language code. See 
##'   \url{https://developers.google.com/+/api/search#available-languages}.
##' @param results The approximate number of results that will be retrieved from
##'   Google+.
##' @param nextToken,cr used internally to retrieve additional pages of answers 
##'   from the Google+ API. Users won't need to set this argument.
##' @return The function returns a list or a data frame containing all available
##'   data on the posts that met the search criteria. See \code{Details} for 
##'   more on its content.
##' @export
##' @seealso Google+ API documentation:
##'   \url{https://developers.google.com/+/api/latest/activities/search}.
##' @examples
##' \dontrun{
##' searchPost("#cats")
##' }
searchPost <- function(q, ret="data.frame", language=NULL, results=1, nextToken=NULL, cr=0) {
  if (is.null(get("apikey", envir=gp))) stop("Set the Google+ API key first using setAPIkey().")
  if (is.null(language)) {
    languageString <- NULL
  } else {
    languageString <- paste0("&language=", language)
  }
  url <- paste0(base.url,
                "activities?query=",
                curlEscape(q),
                languageString,
                "&maxResults=20",
                nextToken,
                "&key=",
                get("apikey", envir=gp))
  this.res <- fromJSON(getURL(url), asText=TRUE)
  res <- this.res[["items"]]
  cr <- cr + length(res)
  if(!is.null(this.res[["nextPageToken"]]) & cr < results) {
    this.nextToken <- paste0("&pageToken=", this.res[["nextPageToken"]])
    res <- c(res, searchPost(q, ret="list", language, results, this.nextToken, cr))
  }
  if (ret=="list") {
    return(res)
  } else {
    res <- ldply(res, parsePost)
    return(res)
  }
}
