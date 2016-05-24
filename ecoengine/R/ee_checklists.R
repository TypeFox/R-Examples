

#' Ecoengine checklists
#'
#' Retrieves exisitng checklists from the ecoengine database
#' @param  subject Enter one of the following subjects: Mammals, Mosses, Beetles, Spiders, Amphibians, Ants, Fungi, Lichen, Plants.
#' @template foptions
#' @export
#' @importFrom assertthat assert_that
#' @importFrom plyr ldply
#' @importFrom httr GET content warn_for_status
#' @return data.frame
#' @examples 
#' # all_lists  <- ee_checklists()
#' # mammals_list  <- ee_checklists(subject = "Mammals")
#' # spiders  <- ee_checklists(subject = "Spiders")
ee_checklists <- function(subject = NULL, foptions = list()) {
	base_url <- paste0(ee_base_url(), "checklists/?format=json")
	full_checklist <- GET(base_url, foptions)
	warn_for_status(full_checklist)
	checklist_data <- content(full_checklist, type = "application/json")
	args <- as.list(ee_compact(c(page_size = checklist_data$count)))
	all_data <- GET(base_url, query = args, foptions)
	warn_for_status(all_data)
	all_checklists <- content(all_data, type = "application/json")
	# all_checklists_df <- rbind_all(all_checklists$results)
	# TO FIX. Remove ldply and see how the speed difference is.
	all_checklists_df <- ldply(all_checklists$results, function(x) data.frame(x))
	if(!is.null(subject)) {
	subject <- eco_capwords(subject)
	sub_result <- all_checklists_df[grep(subject, all_checklists_df$subject), ]
	message(sprintf("Returning %s checklists", nrow(sub_result)))
	return(sub_result)
	} else {
		message(sprintf("Returning %s checklists", all_checklists$count))
		all_checklists_df
	}
}



#'Checklist details
#'
#' Will return details on any checklist 
#' @param list_name URL of a checklist
#' @param  ... Additional arguments (currently not implemented)
#' @export
#' @seealso \code{\link{ee_checklists}}
#' @return \code{data.frame}
#' @examples \dontrun{
#' spiders  <- ee_checklists(subject = "Spiders")
#' # Now retrieve all the details for each species on both lists
#' library(plyr) 
#' spider_details <- ldply(spiders$url, checklist_details)
#'}
checklist_details <- function(list_name, ...) {
details <- GET(paste0(list_name, "?format=json"))
details_data <- content(details, type = "application/json")
first_results <- ldply(details_data$observations, function(x) data.frame(x))
first_results$url <- paste0(first_results$url, "?format=json")
# Now fetch all the results from the URL (2nd column)
full_results <- ldply(first_results$url, function(x) {
				full_checklist <- content(GET(x))
			    rbindfillnull(full_checklist)
})
full_results
}



# Function to capitalize words
#' @noRd
eco_capwords <- function(s, strict = FALSE, onlyfirst = FALSE) {
        cap <- function(s) paste(toupper(substring(s,1,1)),
                {s <- substring(s,2); if(strict) tolower(s) else s}, sep="", collapse=" " )
  if(!onlyfirst){
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
  } else
  {
   sapply(s, function(x) paste(toupper(substring(x,1,1)),tolower(substring(x,2)),sep = "", collapse = " "), USE.NAMES = F)
  }
}
