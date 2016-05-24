
#'Retrieve data from the AntWeb
#'
#' This function allows a user to query the AntWeb database by any taxonomic rank or full species name.
#' @param genus An ant genus name
#' @param  species a species name
#' @param  scientific_name An easier way to pass the Genus and species name together, especially when the data are derived from other packages.
#' @param  georeferenced Default is \code{FALSE}. Set to \code{TRUE} to return only data with lat/long information. Note that this filtering takes place on the client-side, not server side.
#' @param  bbox A lat long bounding box. Format is \code{lat,long,lat,long}. Use this website: http://boundingbox.klokantech.com/ to quickly grab a bbox (set format on bottom left to csv and be sure to switch the order from long, lat, long, lat to lat, long, lat, long)
#' Just set the format on the bottom left to CSV.
#' @param  type A holotype
#' @param  habitat A fuzzy search by any habitat
#' @param  country A country name
#' @param  min_date A lower date bound in the format \code{yyyy-mm-dd}
#' @param  max_date An upper date bound in the format \code{yyyy-mm-dd}
#' @param  min_elevation A lower elevation bound
#' @param  max_elevation An upper elevation bound
#' @param  limit A numeric value to limit number of records
#' @param  offset An offset best used with limit as a way to paginate records
#' @param  quiet If true, any informative messages will be suppressed
#' @export
#' @keywords data download
#' @importFrom rjson fromJSON
#' @importFrom assertthat assert_that
#' @import httr
#' @return data.frame
#' @examples   
#' # data <- aw_data(genus = "acanthognathus", species = "brevicornis")
#' # data3 <- aw_data(genus = "acanthognathus", species = "brevicornis", georeferenced = TRUE)
#' # data2 <- aw_data(scientific_name = "acanthognathus brevicornis")
#' # sandstone <- aw_data(genus = "Aphaenogaster", habitat = "sandstone")
#' # data_genus_only <- aw_data(genus = "acanthognathus", limit = 25)
#' # leaf_cutter_ants  <- aw_data(genus = "acromyrmex")
#' # data  <- aw_data(genus = "Technomyrmex", bbox = '37.77,-122.46,37.76,-122.47')
#' # Search just using a bounding box
#' # data  <- aw_data(bbox = '37.77,-122.46,37.76,-122.47')
#' # Search by a elevation band
#' # aw_data(min_elevation = 1500, max_elevation = 2000)
#' # When you throw a really specimen rich band like below, you'll get a huge number of requests. 
#' # Only the first 1000 records will download first. 
#' # aw_data(min_elevation = 200, max_elevation = 400)
#' # aw_data(min_date = '1980-01-01', max_date = '1981-01-01')
#' # fail <- aw_data(scientific_name = "auberti levithorax") # This should fail gracefully
aw_data <- function(genus = NULL, species = NULL, scientific_name = NULL, georeferenced = NULL, min_elevation = NULL, max_elevation = NULL, type = NULL, habitat = NULL, country = NULL, min_date = NULL, max_date = NULL, bbox = NULL, limit = NULL, offset = NULL, quiet = FALSE) {

	# Check for minimum arguments to run a query
	main_args <- z_compact(as.list(c(scientific_name, genus, type, habitat, bbox)))
	date_args <- z_compact(as.list(c(min_date, max_date)))
	elev_args <- z_compact(as.list(c(min_elevation, max_elevation)))
	arg_lengths <- c(length(main_args), length(date_args), length(elev_args))

	assert_that(any(arg_lengths) > 0)
	decimal_latitude <- NA
	decimal_longitude <- NA
	if(!is.null(scientific_name)) {
		genus <- strsplit(scientific_name, " ")[[1]][1]
		species <- strsplit(scientific_name, " ")[[1]][2]
	}

	# This is a quick call with one result observation requested
	# If all goes well and result set is not greater than 1k, all are retrieved
	base_url <- "http://www.antweb.org/api/v2/"
	original_limit <- limit
	args <- z_compact(as.list(c(genus = genus, species = species, bbox = bbox, min_elevation = min_elevation, max_elevation = max_elevation, habitat = habitat, country = country, type = type, min_date = min_date, max_date = max_date, limit = 1, offset = offset, georeferenced = georeferenced)))
	results <- GET(base_url, query = args)
	warn_for_status(results)
	data <- fromJSON(content(results, "text"))
	data <- z_compact(data) # Remove NULL

	if(data$count > 1000 & is.null(limit)) {
		args$limit <- 1000
		results <- GET(base_url, query = args)
		if(!quiet) message(sprintf("Query contains %s results. First 1000 retrieved. Use the offset argument to retrieve more \n", data$count))
	} else { 
		args$limit <- original_limit
		results <- GET(base_url, query = args)
	}
	
		data <- fromJSON(content(results, "text"))
		data <- z_compact(data)

	if(identical(data$specimens$empty_set, "No records found.")) {
		NULL 
	} else {

	if(!quiet) message(sprintf("%s results available for query.", data$count))
	data_df <- lapply(data$specimens, function(x){ 
	x$images <- NULL	 	
	# In a future fix, I should coerce the image data back to a df and add it here.
	df <- data.frame(t(unlist(x)), stringsAsFactors=FALSE)
	df
})
	final_df <- data.frame(do.call(rbind.fill, data_df))
	names(final_df)[grep("geojson.coord1", names(final_df))] <- "decimal_latitude"
	names(final_df)[grep("geojson.coord2", names(final_df))] <- "decimal_longitude"
	# There seem to be extra field when searching for just a genus
	final_df$decimalLatitude <- NULL
	final_df$decimalLongitude <- NULL
	final_df$minimumElevationInMeters <- as.numeric(final_df$minimumElevationInMeters)
	final_results <- list(count = data$count, call = args, data = final_df)
	class(final_results) <- "antweb"
	final_results

}
}	
 

#' Download all aw_data available for any request
#'
#' This is a thin wrapper around aw_data
#' @param ... All the same arguments that get passed to \code{aw_data}
#' @param progress Default is on and set to \code{text}. Set to \code{none} to suppress
#' @export
#' @importFrom plyr llply
#' @keywords data download
#' @seealso aw_data
#' @examples \dontrun{
#' # crem <- aw_data_all(genus = "crematogaster", georeferenced = TRUE)
#'}
aw_data_all <- function(..., progress = 'text') {
	x <- aw_data(..., quiet = TRUE)
	message(sprintf("Downloading %s results. Be patient\n", x$count))
	if(!is.null(x$count)) {
	bins <- seq(from = 0, to = x$count, by = 1000)
	results <- llply(bins, function(x) {
		aw_data(..., offset = x, quiet = TRUE)
	}, .progress = progress)

	aw_cbind(results)
}
} 



#' aw_unique
#'
#' Get a list of unique names within any taxonomic rank
#' @param rank  A taxonomic rank. Allowed values are  \code{subfamily}, \code{genus} or \code{species}
#' @param  name Optional. If left blank, the query will return a list of all unique names inside the supplied rank.
#' @export
#' @seealso \code{\link{aw_data}}
#' @return data.frame
#' @examples  \dontrun{
#' subfamily_list <- aw_unique(rank = "subfamily")
#' # genus_list <- aw_unique(rank = "genus")
#' # species_list <- aw_unique(rank = "species")
#'}
aw_unique <- function(rank = NULL, name = NULL) {
	# assert_that(!is.null(z_compact(c(rank, name))))
	# base_url <- "http://www.antweb.org/api/v2"
	# args <- z_compact(as.list(c(rank = rank, name = name)))
	# results <- GET(base_url, query = args)
	# warn_for_status(results)
	# data <- fromJSON(content(results, "text"))
	# data.frame(do.call(rbind, data))
	message("Function deprecated. Use aw_distinct instead. Or maybe not. Might fold in aw_distnct here.")
}

