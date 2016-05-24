#' Calculates the area of a given region.
#' 
#' Calculates the area of a given region to a given precision.
#' 
#' 
#' @param data The region, should contain \$lat and \$lon.
#' @param ngrdpts The precision of the calculation.
#' @param Projection Which projection is beeing used.
#' @param old.method If true an older version of this program is used. Default
#' is False.
#' @param robust If true a more robust method is used, default is True.
#' @return The area of the region given.
#' @section Side Effects: None
#' @seealso \code{\link{geodefine}}, \code{\link{geolocator}},
#' \code{\link{geoinside}}.
#' @examples
#' 
#'          geoarea(island)         # Calculates the area of Iceland up to 
#'                                  # an with default precision.
#' 
#'          geoarea(island,10000)   # Calculates the area of Iceland up to 
#'                                  # an adiquite precision.
#' 
#' #         geoarea(geodefine(),10) # Calculates the area of a region specified 
#'                                  # by the user.
#' 
#' @export geoarea
geoarea <-
function(data, Projection = "Lambert", old.method = F, ngrdpts = 2000, robust
	 = T)
{
	area <- 0
	data <- geo.Split.poly(data)
	if(old.method) {
		for(i in 1:length(data))
			area <- area + geoarea.old(data[[i]], ngrdpts, robust)
	}
	else {
		area <- 0
		for(i in 1:length(data)) {
			if(Projection == "Lambert")
				data[[i]] <- lambert(data[[i]]$lat, data[[i]]$
					lon, mean(data[[i]]$lat), mean(data[[
					i]]$lon), mean(data[[i]]$lat))
			else data[[i]] <- mercator(data[[i]]$lat, data[[i]]$
					lon, b0 = mean(data[[i]]$lat))
			data[[i]] <- data.frame(x = data[[i]]$x, y = data[[
				i]]$y)
			n <- nrow(data[[i]])
			area <- area + abs(sum(data[[i]]$x[1:(n - 1)] * data[[
				i]]$y[2:n] - data[[i]]$x[2:n] * data[[i]]$
				y[1:(n - 1)], na.rm = T)/2)
		}
	}
	return(area)
}

