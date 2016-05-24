#' An tematres api helper that converts ids to strings and vise versa
#'
#' As the tematres api accepts ids of terms only on most of its tasks this is a small
#' helper that converts them into each other. It has two function names tha you can
#' use depending on what you are converting. But both function names do the same so
#' they can be exchanged. (rtematres.api.conversion.[term_id/id_term])
#'
#' @param given Either a term or an integer id for a term
#' @param warn Give a waring when the term/id is not found. Defaults to TRUE.
#'
#' @return Either a string or id
#' @export

rtematres.api.conversion.term_id <- rtematres.api.conversion.id_term <- function(given, warn = T) {
	if(is.character(given)) {
		response = rtematres.api(task = "search", argument = given)
		the_id_is = response$id[which(response$term == given)]
		if(length(the_id_is) == 0 && warn) warning("Sorry no such term available!")
		return(the_id_is)
	} else {
		response = rtematres.api(task = "fetchTerm", argument = given)
		the_term_is = response$term[which(response$id == given)]
		if(length(the_term_is) == 0 && warn) warning("Sorry no such id available!")
		return(the_term_is)
	}
}
