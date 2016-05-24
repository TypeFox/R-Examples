#' Create spatial data object for model fitting via \code{spatialOccupancy}
#' function
#' 
#' This function takes an observation data frame and a data frame of site
#' characteristics and combines them together for analysis with the
#' \code{\link{spatial.occupancy}} function.
#' 
#' This function combines the two data frames and assigns names so that
#' \code{spatial.occupancy} knows which columns to use. It also performs
#' some rudimentary error checking to make sure the data is in the proper form
#' (e.g., the site IDs in the visit data frame must be contained in the site
#' IDs for the site data frame)
#' 
#' @param visit.data A data frame that contains the observed occupancy for each
#' site and any detection related covariates.
#' @param site.data A data frame that contains the site id, coordinates, and
#' any habitat related covariates that might influence the occupancy process
#' @param names A named list with the following elements: (1)\code{visit} A
#' named list with elements "site" = the name of the site id in the observation
#' data frame and "obs" = the name of the observed occupancy variable (2)
#' \code{site} A named list with elements "site" = the name of the site id and
#' "coords" = a character vector giving the name of the coordinates (x first
#' then y)
#' @return An \code{so.data} object is a list with elements equal to the two
#' data frames. Attributes are set giving the names of columns of interest
#' @author Devin S. Johnson <devin.johnson@@noaa.gov>
#' @export
make.so.data <-
function(visit.data, site.data, names){
	visit.names <- names$visit
	site.names <- names$site
	if(is.null(visit.data[,visit.names$site])) stop(visit.names$site," is not in the visit data frame!")
	if(is.null(site.data[,site.names$site])) stop(site.names$site," is not in the site data frame!")
	if(is.null(visit.data[,visit.names$obs])) stop(visit.names$obs," is not in the visit data frame!")
	if(is.null(visit.data[,visit.names$site])) stop(visit.names$site," is not in the visit data frame!")
	if(any(!visit.data[,visit.names$obs]%in%c(0,1))) stop("There are observed detection values other than 0 or 1!")
	if(any(!visit.data[,visit.names$site]%in%site.data[,site.names$site])) stop("There are visited sites not listed in the site data frame!")
	attr(visit.data, "site") <- visit.names$site
	attr(visit.data, "obs") <- visit.names$obs
	attr(site.data, "site") <- site.names$site
	attr(site.data, "coords") <- site.names$coords
	so.data <- list(visit=visit.data, site=site.data)
	class(so.data) <- "so.data"
	return(so.data)
}

