
#' LeafletJS Map
#'
#' Builds an interactive map of locations for any list of species
#' @param aw_obj Result from a search on AntWeb
#' @param  dest Location where the html file and geojson file should be stored. Default is the temp directory
#' @param  title Title of the map.
#' @param  incl.data Default is \code{TRUE}. Writes geoJSON data into the html file to get around security restrictions in browsers like Google Chrome. Set to \code{FALSE} to read from a separate local geoJSON file.
#' @export
#' @keywords map
#' @import leafletR
#' @examples \dontrun{
#'  acanthognathus_df <- aw_data(genus = "acanthognathus", georeferenced = TRUE)
#'  aw_map(acanthognathus_df)
#' # Or just plot data by habitat. So for e.g. using sandstone as a substrate
#' sandstone <- aw_data(habitat = "sandstone")
#' aw_map(sandstone)
#'}
aw_map <- function(aw_obj, dest = tempdir(), title = "AntWeb species map", incl.data = TRUE) {
	decimal_latitude <- NULL
	decimal_longitude <- NULL

	assert_that(identical(class(aw_obj), "antweb"))
	aw_obj <- aw_obj$data
	aw_obj <- subset(aw_obj, !is.na(decimal_latitude) & !is.na(decimal_longitude))
	assert_that(nrow(aw_obj) >= 1)


	dest <- ifelse(is.null(dest), tempdir(), dest)
	aw_obj$scientific_name <- aw_obj$scientific_name
	species_data <- aw_obj
	lat_location <- which(names(species_data) == "decimal_latitude")
	lon_location <- which(names(species_data) == "decimal_longitude")
	num_species <- length(unique(species_data$scientific_name))
	# 	if(num_species > 30) {
	# 		stop("Map cannot accommodate more than 30 species. Please choose a smaller subset of data.")
	# }

	ee_geo <- toGeoJSON(data = species_data, name = "temp", dest = dest, lat.lon = c(lat_location, lon_location))	
	if(length(unique(species_data$scientific_name)) <= 30) {
	cols <- c("#8D5C00", "#2F5CD7","#E91974", "#3CB619","#7EAFCC",
"#4F2755","#F5450E","#264C44","#3EA262","#FA43C9","#6E8604","#631D0E","#EE8099","#E5B25A",
"#0C3D8A","#9E4CD3","#195C7B","#9F8450","#7A0666","#BBA3C5","#F064B4","#108223","#553502",
"#17ADE7","#83C445","#C52817","#626302","#9F9215","#6CCD78","#BF3704")
	pal <- cols[1:num_species]
	sty <- styleCat(prop = "scientific_name", val = unique(species_data$scientific_name), style.val = pal, fill.alpha = 1, alpha = 1, rad = 4, leg = "Scientific Name")
	map <- leaflet(ee_geo, base.map="tls", style = sty, popup = "scientific_name", dest = dest, title = title, incl.data = incl.data) 
	} else {
	message("\nDataset contains more than 30 unique species so a legend was not generated\n")	
	map <- leaflet(ee_geo, base.map="tls", dest = dest, title = title, incl.data = incl.data) 		
}
	browseURL(map)
}

