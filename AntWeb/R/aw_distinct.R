
#' aw_distinct
#'
#' Retrieves a data.frame of distinct ranks based on various restrictions
#' @param rank = "genus" Default is genus. But you can also use phylum, sub-phylum etc
#' @param  habitat  The habitat type
#' @param  country Country name
#' @param  min_elevation  Min elevation recorded for specimen
#' @param  max_elevation  Max elevation recorded for specimen
#' @param  limit = 1000 Default limit. Set higher if necessary
#' @param  offset To be used in conjunction with limit
#' @export
#' @examples \dontrun{
#' aw_distinct(rank = "genus", country = "Madagascar")
#'}
aw_distinct <- function(rank = "genus", 
						habitat = NULL, 
						country = NULL, 
						min_elevation = NULL, 
						max_elevation = NULL, 
						limit = 1000, 
						offset = NULL) {
	main_args <- z_compact(as.list(c(rank = rank, 
									habitat = habitat, 
									country = country, 
									min_elevation = min_elevation, 
									max_elevation =  max_elevation, 
									limit = limit, 
									offset = offset)))
	base_url <- aw_base_url()
	results <- GET(base_url, query = main_args)
	warn_for_status(results)
	data <- fromJSON(content(results, "text"))
	results <- data.frame(unlist(data))
	names(results) <- rank
	final_results <- list(count = data$count, call = main_args, data = results)
	class(final_results) <- "antweb"
	final_results

}