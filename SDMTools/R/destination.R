#' Vincenty Direct Calculation of a Destination
#' 
#' \code{destination} estimates the destination latitude and longitude given a
#' starting latitude and longitude, a bearing and distance. \cr \cr For general
#' information on Vincenty's formula, see e.g.,
#' \url{http://en.wikipedia.org/wiki/Vincenty's_formulae}. It states: \cr
#' \emph{Vincenty's formulae are two related iterative methods used in geodesy
#' to calculate the distance between two points on the surface of an spheroid,
#' developed by Thaddeus Vincenty in 1975. They are based on the assumption
#' that the figure of the Earth is an oblate spheroid, and hence are more
#' accurate than methods such as great-circle distance which assume a spherical
#' Earth.} \cr \cr \bold{Note:} this method assumes a locations are lat & lon
#' given in WGS 84.
#' 
#' Typical useages are:\cr \enumerate{ \item a single start location, bearing
#' and distance to give a single output location\cr --output would be a single
#' destination location \item a single start location with one or more bearings
#' or distances to give multiple output locations\cr --output would be a
#' destination locations for each combination of bearings and distances \item
#' multiple start locations with a single bearing or distance\cr --output would
#' be a destination locations representing the bearing and distance from each
#' of the start locations \item multiple start locations with multiple bearings
#' or distances\cr --output would be a destination locations representing the
#' combinations of bearings and distances from each of the start locations\cr
#' -- NOTE that the bearing and distance vectors must be of the same length of
#' the input lat and long.  }
#' 
#' See examples for all possible usages.
#' 
#' @param lat a single value or vector of values representing latitude in
#' decimal degrees from -90 to 90 degrees.
#' @param lon a single value or vector of values representing longitude in
#' decimal degrees from -180 to 180 degrees.
#' @param bearing a single value or vector of values representing the bearings
#' (directions) of interest ranging from 0 to 360 degrees.
#' @param distance a single value or vector of values representing the
#' distances in metres to the destination.
#' @return Returns a data.frame with: \item{lon1}{the original longitude}
#' \item{lat1}{the original latitude} \item{bearing}{the bearing used}
#' \item{distance}{the distance used} \item{lon2}{the destination longitude}
#' \item{lat2}{the destination latitude}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @references Vincenty, T. 1975. Direct and Inverse Solutions of Geodesics on
#' the Ellipsoid with application of Nested Equations. Survey Review, vol XXII
#' no 176. \url{http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf}
#' @source The source code here was modified from
#' \url{http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html}.\cr
#' \cr Destinations were validated against Geoscience Australia calculations
#' (\url{http://www.ga.gov.au/geodesy/datums/vincenty_direct.jsp}).
#' @examples
#' 
#' 
#' ###single lat lons
#' lats = -85; lons = 165
#' #single bearing & single distance
#' destination(lats,lons,bearing=180,distance=500000)
#' 
#' #multiple bearings
#' destination(lats,lons,bearing=seq(0,360,length.out=9),distance=500000)
#' 
#' #multiple bearings
#' destination(lats,lons,bearing=45,distance=seq(0,5000000,length.out=11))
#' 
#' #multiple bearings, multiple distances
#' destination(lats,lons,bearing=seq(0,360,length.out=9),
#'     distance=seq(0,5000000,length.out=11))
#' 
#' ###multiple lat lons
#' lats = seq(-90,90,length.out=9); lons = seq(-180,180,length.out=9)
#' 
#' #multiple lat lons but single bearings / distances
#' destination(lats,lons,bearing=45,distance=500000)
#' 
#' #different bearings for each lat lon
#' destination(lats,lons,bearing=seq(0,360,length.out=9),distance=500000)
#' 
#' #different distances for each lat lon
#' destination(lats,lons,bearing=45,distance=seq(0,5000000,length.out=9))
#' 
#' #different bearings & distances for each lat lon
#' destination(lats,lons,bearing=seq(0,360,length.out=9),
#'     distance=seq(0,5000000,length.out=9))
#' 	
#' 
#' @export 
#' @useDynLib SDMTools Dest
destination = function(lat, lon, bearing, distance) {
	#check the data
	if (length(lat)!=length(lon)) stop('lat & lon must be of the same length')
	if (any(lon < -180) | any(lon > 180)) stop('lon must be decimal degrees between 0 & 360')
	if (any(lat < -90) | any(lat > 90)) stop('lat must be decimal degrees between -90 & 90')
	if (length(lat)==1) {
		out = expand.grid(bearing=bearing,distance=distance)
		out = data.frame(lon1=lon,lat1=lat,out,lon2=NA,lat2=NA)
	} else {
		if (length(bearing)>1) { if (length(bearing)!=length(lat)) stop('number of bearing values is not the same length as lon & lat') } else { bearing = rep(bearing,length(lat)) }
		if (length(distance)>1) { if (length(distance)!=length(lat)) stop('number of distance values is not the same length as lon & lat') } else { distance = rep(distance,length(lat)) }
		out = data.frame(lon1=lon,lat1=lat,bearing=bearing,distance=distance,lon2=NA,lat2=NA)
	}
	#clatcle through and output the new data
	for (ii in 1:nrow(out)) {
		tt = .Call('Dest',out$lat1[ii],out$lon1[ii],out$bearing[ii],out$distance[ii],PACKAGE='SDMTools')
		out$lon2[ii] = tt[2]; out$lat2[ii] = tt[1]
	}
	#return the output
	return(out)
}
