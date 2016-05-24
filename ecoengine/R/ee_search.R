

#' Ecoengine search
#'
#' Search across the entire ecoengine database.
#' @param query search term
#' @template foptions
#' @importFrom dplyr rbind_all
#' @export
#' @keywords search
#' @examples \dontrun{
#' # lynx_results <- ee_search(query = "genus:Lynx")
#'}
ee_search <- function(query = NULL, foptions = list()) {

# search_url <- "http://ecoengine.berkeley.edu/api/search/?format=json"
search_url <- paste0(ee_base_url(), "search/?format=json")
args <- as.list(compact(c(q = query)))
result <- GET(search_url, query = args, foptions)
es_results <- content(result, type = "application/json")
fields <- es_results$fields
# This removes list items with nothing nested.
ee_filter <- function(i) {
        length(i) > 0
}
fields_compacted <- Filter(ee_filter, fields)
faceted_search_results <- lapply(fields_compacted, function(y) {
        temp_fields <- do.call(rbind.data.frame, lapply(y, LinearizeNestedList))
        names(temp_fields) <- c("field", "results", "search_url")
        temp_fields
})
do.call(rbind, faceted_search_results)
}



#'Search observations
#'
#' A powerful way to search through the observations.
#' @param query = The search term
#' @template foptions
#' @param  quiet Default is \code{FALSE}. Set to \code{TRUE} to suppress messages.
#' @template pages
#' @template progress
#' @export
#' @keywords search
#' @seealso \code{\link{ee_search})}
#' @return data.frame
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' # general_lynx_query <- ee_search_obs(query  = "Lynx")
#'  \dontrun{
#' # lynx_data <- ee_search_obs(query  = "genus:Lynx")
#' # all_lynx_data <- ee_search_obs(query  = "Lynx", page = "all")
#'}
ee_search_obs <- function(query = NULL, page = NULL, page_size = 1000, quiet = FALSE, progress = TRUE, foptions = list()) {
	# obs_search_url <- "http://ecoengine.berkeley.edu/api/observations/?format=json"
	obs_search_url <- paste0(ee_base_url(), "observations/?format=json")
	args <- ee_compact(as.list(c(q = query, page_size = 1000)))
    main_args <- args
    main_args$page <- page
	obs_search <- GET(obs_search_url, query = args, foptions)
	obs_results <- content(obs_search, type = "application/json")
	if(is.null(page)) { page <- 1 }
	required_pages <- ee_paginator(page, obs_results$count)

if(!quiet) {
message(sprintf("Search contains %s observations (downloading %s of %s pages)", obs_results$count, length(required_pages), max(required_pages)))
}

if(progress) pb <- txtProgressBar(min = 0, max = length(required_pages), style = 3)
    results <- list()
    for(i in required_pages) {
    	args$page <- i
    		obs_search <- GET(obs_search_url, query = args, foptions)
			obs_results <- content(obs_search, type = "application/json")
			obs_results2 <- lapply(obs_results$results, LinearizeNestedList)
			obs_data <- lapply(obs_results2, function(x) data.frame(t(unlist(x))))
			results[[i]] <- obs_data
			if(progress) setTxtProgressBar(pb, i)
		}
	all_obs_data <- do.call(rbind.fill, results)
	all_obs_results <- list(results = obs_results$count, call = main_args, type = "observations", data = all_obs_data)
	class(all_obs_results) <- "ecoengine"

	if(progress) close(pb)
    all_obs_results
}





# @importFrom data.table rbindlist
