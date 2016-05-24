#' Advanced Search.  
#'
#' Search function that will filter on matching timestamp, author, title, description, tag, category, and date range.  Query searches against matches in any metadata field.  Full-text searches coming soon.  
#' @author Carl Boettiger \email{cboettig@@gmail.com}
#' @param query the search query
#' @param author Show only results by this author
#' @param title Show only results matching or partially matching this title
#' @param description Show only results matching or partially matching this description
#' @param tag Show only results matching this tag
#' @param category Show only results matching this category
#' @param from_date Start time window for search. Date format is YYYY-MM-DD
#' @param to_date Ending time window for search. Date format is YYYY-MM-DD
#' @param mine Browse only articles owned by user. default is FALSE. Not functional. Use \code{\link{fs_browse}} instead.   
#' @param public_only (for use with mine=TRUE only) browse only my public articles. default is FALSE
#' @param private_only (for use with mine=TRUE only) browse only my private articles. default is FALSE
#' @param drafts_only (for use with mine=TRUE only) browse only my draft articles. default is FALSE
#' @param session (optional) the authentication credentials from \code{\link{fs_auth}}. If not provided, will attempt to load from cache as long as figshare_auth has been run.  
#' @param base the API access url
#' @param debug logical, enable debugging mode
#' @return output of PUT request (invisibly)
#' @seealso \\code{\link{fs_auth}} \code{\link{fs_browse}}
#' @references \url{http://api.figshare.com/docs/howto.html#q-search}
#' @import httr
#' @export 
#' @examples \dontrun{
#' fs_search("Boettiger") 
#' fs_search("Boettiger", author = "Carl")
#' fs_search("Boettiger", author = "Carl", from="2014-01-01")
#' fs_search("Boettiger", author = "Carl", from="2014-01-01",
#'           category = "Evolutionary Biology")
#' 
#' } 
fs_search <- 
  function(query, 
           author = NA, 
           title = NA, 
           description = NA, 
           tag = NA, 
           category = NA, 
           from_date = NA, 
           to_date = NA, 
           mine = FALSE, 
           public_only = FALSE, 
           private_only = FALSE, 
           drafts_only = FALSE, 
           session = fs_get_auth(), 
           base = "http://api.figshare.com/v1", 
           debug = FALSE) {

    # resolve name conflict?
    the_title <- title

    ## Fine assembling queries this way, but would be prettier
    ## using parameters to RCurl as in RMendeley...
    method <- "articles"
    if(mine)
      method <- "my_data/articles"
    if(public_only)
      method <- paste(method, "/public", sep="")  # visibility only works in my_data
    if(private_only)
      method <- paste(method, "/private", sep="") # visibility only works in my_data
    if(drafts_only)
      method <- paste(method, "/drafts", sep="")  # visibility only works in my_data
    if(!is.na(query)) # Can search work in my_data ?  
      method <- paste(method, "/search?search_for=", query, sep="")


    if(!is.na(author)){ # NOT AN ID
      method <- paste(method, "&has_author=", author, sep="")
    }
    if(!is.na(the_title))
      method <- paste(method, "&has_title=", the_title, sep="")
    if(!is.na(description))
      method <- paste(method, "&has_description=", description, sep="")
    if(!is.na(tag))
      method <- paste(method, "&has_tag=", tag, sep="")
    if(!is.na(category)){
#      if(!is.numeric(category))
#        category <- fs_cat_to_id(category)
      method <- paste(method, "&has_category=", category, sep="")
    } 
    if(!is.na(from_date))
      method <- paste(method, "&from_date=", from_date, sep="")
    if(!is.na(to_date))
      method <- paste(method, "&to_date=", to_date, sep="")

    request <- paste(base, method, sep="/")
    request <- build_url(parse_url(request)) # perform % encoding

    out <- GET(request, config(token = session))
    if(debug | out$status_code != 200)
      out

    else {

      parsed <- fromJSON(content(out, "text"))
      if(is.null(parsed$count))
        parsed$count <- parsed$items_found
      if(is.null(parsed$count))
        parsed$count <- length(parsed$items)
      out <- parsed$items



      if(parsed$count > 10) {
        total_pages <- ceiling(parsed$count / 10)
        all <- lapply(1:total_pages, function(i){
          method_ <- paste(method, "&page=", i, sep="")
          request = paste(base, method_, sep="/")
          request <- build_url(parse_url(request)) # perform % encoding
          out <- GET(request, config(token = session))
          parsed <- fromJSON(content(out, as = "text"))
          parsed$items
        })
        out <- unlist(all, recursive = FALSE)
      }
#      out <- lapply(out, function(o){
#             class(o) <- "fs_object"
#             o 
#        })
      out
    }
  }

## response <- GET("http://api.figshare.com/v1/articles/search?search_for=Boettiger", config(token = fs_auth())) 


