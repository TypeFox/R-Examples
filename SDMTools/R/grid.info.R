#' Grid Information from Geographic (lat lon) Projections
#' 
#' Since spatial grids in geographic projections do not have equal area or
#' perimeters, \code{grid.info} extracts perimeter & area related information
#' for latitudinal bands with differing longitudinal widths. \cr\cr Outputs
#' lengths are in m using Vincenty's equation (\code{distance})and areas in m2.
#' Surface areas are calculated summing surface areas of spherical polygons as
#' estimated using l'Huiller's formula.
#' 
#' 
#' @param lats is a vector of latitudes representing the midpoint of grid cells
#' @param cellsize is a single value (assuming square cells) or a two value
#' vector (rectangular cells) representing the height (latitude) and width
#' (longitude) of the cells
#' @param r is a single value representing the radius of the globe in m.
#' Default is for the WGS84 elipsoid
#' @return a data.frame listing: \item{lat}{the latitude representing the
#' midpoint of the cell} \item{top}{length of the top of the cell (m)}
#' \item{bottom}{length of the bottom of the cell (m)} \item{side}{length of
#' the side of the cell (m)} \item{diagnal}{length of the diagnals of the cell
#' (m)} \item{area}{area of the cell (m2)}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @references information on l'Huiller's formula
#' \url{http://williams.best.vwh.net/avform.htm for more info)} code for
#' estimating area of polygon on sphere was modified from
#' \url{http://forum.worldwindcentral.com/showthread.php?t=20724}
#' @examples
#' 
#' #show output for latitudes from -87.5 to 87.5 at 5 degree intervals
#' grid.info(lats=seq(-87.5,87.5,5), 5)
#' 
#' @export 
grid.info = function(lats,cellsize,r=6378137) {
	r2 = r^2 #radius of earth
	###need checks to ensure lats will not go beyond 90 & -90
	if (length(cellsize)==1) cellsize=rep(cellsize,2) #ensure cellsize is defined for both lat & lon
	out = data.frame(lat=lats) #setup the output dataframe
	toplats = lats+(0.5*cellsize[1]); bottomlats = lats-(0.5*cellsize[1]) #define the top and bottom lats
	check = range(c(toplats,bottomlats),na.rm=TRUE); if (-90>check[1] | 90<check[2]) stop('latitudes must be between -90 & 90')
	out$top = distance(toplats,rep(0,length(lats)),toplats,rep(cellsize[2],length(lats)))$distance
	out$bottom = distance(bottomlats,rep(0,length(lats)),bottomlats,rep(cellsize[2],length(lats)))$distance
	out$side = distance(toplats,rep(0,length(lats)),bottomlats,rep(0,length(lats)))$distance
	out$diagnal = distance(toplats,rep(0,length(lats)),bottomlats,rep(cellsize[2],length(lats)))$distance
	#calculate area of a spherical triangle using spherical excess associated by knowing distances
	#tan(E/4) = sqrt(tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2))
	#where a, b, c = sides of spherical triangle
	#s = (a + b + c)/2
	#from CRC Standard Mathematical Tables
	#calculate excess based on  l'Huiller's formula (http://williams.best.vwh.net/avform.htm for more info)
	#code modified from (http://forum.worldwindcentral.com/showthread.php?t=20724)
	excess = function(lam1,lam2,beta1,beta2){ #calculate excess... inputs are in radians
		haversine = function(y) { (1-cos(y))/2 }
		cosB1 = cos(beta1); cosB2 = cos(beta2)
		hav1 = haversine(beta2-beta1) + cosB1*cosB2*haversine(lam2-lam1)
		aa = 2 * asin(sqrt(hav1)); bb = 0.5*pi - beta2; cc = 0.5*pi - beta1
		ss = 0.5*(aa+bb+cc)
		tt = tan(ss/2)*tan((ss-aa)/2)*tan((ss-bb)/2)*tan((ss-cc)/2)
		return(abs(4*atan(sqrt(abs(tt)))))		
	}
	if (any(bottomlats==-90)) { pos = which(bottomlats==-90); bottomlats[pos] = -bottomlats[pos]; toplats[pos] = -toplats[pos]} #ensure no -90 bottom lats
	out$area = excess(lam1=0,lam2=cellsize[2]*pi/180,toplats*pi/180,toplats*pi/180)
	out$area = abs(out$area-excess(lam1=0,lam2=cellsize[2]*pi/180,bottomlats*pi/180,bottomlats*pi/180))*r2
	return(out)
}
