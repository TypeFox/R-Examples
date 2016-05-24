##' Retrieve the posts of a user's G+ page
##' 
##' This function retrieves the most recent posts that a user put on her 
##' page. Google calls this `listing activities`.
##' 
##' The result is either a simple list of items from the page that can be parsed
##' using e.g. \code{\link{parsePost}} or a data frame with that (or another 
##' user-supplied) function already applied.
##' 
##' When writing your own parsing functions, make sure that the function takes a
##' single list item from the raw list as its argument and returns a vector of
##' values or a one-row data frame. The return values of the function are then
##' fed into \code{plyr}'s \code{ldply} to turn it into a data frame. See
##' \code{\link{parsePost}} for an example.
##' 
##' The length of the list or the number of rows of the data frame are somewhat 
##' ambiguous. Specifying the \code{results} argument will try to get that many
##' results. But there may be less (because Google could not find more) or more 
##' (because Google is organizing results on pages and it would be a waste to 
##' discard them automatically). If you really depend on getting not more rows 
##' than you expected, use standard selection (i.e. \code{[}) to trim the
##' results.
##' 
##' @param user A user identification string: either user ID or +Name.
##' @param parseFun A function for parsing the results, e.g. the supplied 
##'   \code{\link{parsePost}} function. If set to \code{NULL}, then the raw list
##'   of the retrieved posts is being returned. Defaults to \code{parsePost}.
##' @param results The approximate number of results that will be retrieved from
##'   Google+.
##' @param nextToken,cr used internally to retrieve additional pages of
##'   answers from the Google+ API. Users won't need to set these arguments.
##' @return The function returns either the raw list of retrieved posts or
##'   whatever the supplied parsing function does with the retrieved list.
##' @export
##' @seealso Google+ API documentation:
##'   \url{https://developers.google.com/+/api/latest/activities/list}.
#' @importFrom plyr ldply
##' @examples
##' \dontrun{
##' myPosts.df <- harvestPage("115046504166916768425")
##' gPosts.df <- harvestPage("+google", results=200)
##' }
harvestPage <- function(user, parseFun=parsePost, results=1, nextToken=NULL, cr=1) {
  if (is.null(get("apikey", envir=gp))) stop("Set the Google+ API key first using setAPIkey().")
  if (results < 1) stop("Argument 'results' needs be positive.")
  url <- paste0(base.url,
                start.people,
                curlEscape(user),
                close.page1,
                nextToken,
                close.page2,
                get("apikey", envir=gp))
  this.res <- fromJSON(getURL(url), asText=TRUE)
  res <- this.res[["items"]]
  cr <- cr + length(res)
  if(!is.null(this.res[["nextPageToken"]]) & cr < results) {
    this.nextToken <- paste0("&pageToken=", this.res[["nextPageToken"]])
    res <- c(res, harvestPage(user, NULL, results, this.nextToken, cr))
  }
  if (is.null(parseFun)) {
    return(res)
  } else {
    res <- ldply(res, parseFun)
    return(res)
  }
}


##' Retrieve the users that acted on a G+ post
##'
##' This function retrieves the users that either +1ed or reshared a post.
##' Google calls this `list by activity`.
##'
##' @param activity The post ID for which the users should be retrieved.
##' @param kind Denoting the kind of person to be retrieved. Either
##'   \code{plusoners} or \code{resharers}.
##' @param nextToken This is used internally to retrieve additional pages of
##'   answers from the Google+ API. Users won't need to set this argument.
##' @return Returns a (character) vector of Google+ user IDs.
##' @seealso Google+ API documentation:
##'   \url{https://developers.google.com/+/api/latest/people/listByActivity}.
##' @export
##' @examples
##' \dontrun{
##' ## User IDs of people that +1ed this post
##' users.p <- harvestActivity("z131ihvycs30ivrxm04cjbiwjkbqujka0sk0k", "plusoners")
##' }
harvestActivity <- function(activity, kind=c("plusoners", "resharers"),
                            nextToken=NULL) {
  if (is.null(get("apikey", envir=gp))) stop("Set the Google+ API key first using setAPIkey().")
  if (kind != "plusoners" & kind != "resharers") stop("Argument 'kind' needs to be either 'plusoners' or 'resharers'.")
  this.url <- paste0(base.url, "activities/",
                     curlEscape(activity),
                     "/",
                     start.people,
                     kind,
                     "?maxResults=100",
                     nextToken,
                     "&key=",
                     get("apikey", envir=gp))
  this.res <- fromJSON(getURL(this.url), asText=TRUE)
  this.ppl <- sapply(this.res[["items"]], function(x) x$id)
  if (!is.null(this.res[["nextPageToken"]])) {
    this.nextToken <- paste0("&pageToken=", this.res[["nextPageToken"]])
    this.ppl <- c(this.ppl,
                  harvestActivity(activity, kind, this.nextToken))
  }
  return(this.ppl)
}


##' Retrieve the profile of Google+ users
##'
##' This function retrieves the profile of one or more Google+ user(s). Google
##' calls this `get people`. The results are either returned as a raw list with 
##' one element per profile or parsed using a parsing function, either the 
##' prepackaged one \code{\link{parseProfile}} or a user-supplied one.
##' 
##' When using your own parsing function, be sure that it takes a single element
##' from the returned list and returns either a vector of values or a single row
##' data frame.
##'
##' @param id A character vector of the Google+ user ID(s).
##' @param parseFun the function used to parse the results. If \code{NULL} the
##' raw list of results is returned.
##' @return The function returns either a raw list or a parsed version. 
##' See \code{Details}.
##' @seealso Google+ API documentation:
##'   \url{https://developers.google.com/+/api/latest/people/get}
##' @export
##' @examples
##' \dontrun{
##' gProfile <- harvestProfile("+google")
##' }
harvestProfile <- function(id, parseFun=parseProfile) {
  if (is.null(get("apikey", envir=gp))) stop("Set the Google+ API key first using setAPIkey().")
  this.res <- lapply(as.list(id), function(x) {
    u <- paste0(base.url,
               start.people,
               curlEscape(x),
               close.people,
               get("apikey", envir=gp))
    r <- fromJSON(getURL(u), asText=TRUE)
    return(r)
  })
  if (is.null(parseFun)) {
    return(this.res)
  } else {
    return(ldply(this.res, parseFun))  
  }
}
