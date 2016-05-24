#' aw_code
#'
#' Retrieve data by specimen id
#' @param occurrenceid A unique id in the AntWeb database identifying a particular specimen
#' @param catalogNumber Specimen catalogue number
#' @export
#' @seealso \code{\link{aw_data}}
#' @return list
#' @examples 
#' # data_by_code <- aw_code(occurrenceid = "CAS:ANTWEB:alas188691") 
#' # data_by_code <- aw_code(catalognumber="inb0003695883")

aw_code <- function(occurrenceid = NULL, catalogNumber = NULL) {

	# We need at least one identifier
	assert_that(!is.null(occurrenceid) | !is.null(catalogNumber))

	occurrenceid <- tolower(occurrenceid)
	base_url <- "http://www.antweb.org/api/v2"
	args <- z_compact(as.list(c(occurrenceId = occurrenceid, catalogNumber = catalogNumber)))
	results <- GET(base_url, query = args)
	warn_for_status(results)
	data <- fromJSON(content(results, "text"))
	if(identical(data$specimens$empty_set, "No records found.")) {
		NULL 
	} else {
	final_df <- data.frame(t(unlist(data[2])))
	names(final_df)[grep("geojson.coord1", names(final_df))] <- "decimal_latitude"
	names(final_df)[grep("geojson.coord2", names(final_df))] <- "decimal_longitude"
	final_results <- list(count = data$count, call = args, data = final_df)
	class(final_results) <- "antweb"
	final_results
}
}
 