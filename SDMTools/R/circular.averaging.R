#' Circular Averaging based on Vector Averaging
#' 
#' \code{circular.averaging} calculates the average direction (0 - 360) given a
#' vector of directions.\cr\cr \code{vector.averaging} calculates the average
#' distance and direction given a vector of directions and a vector of
#' distances.
#' 
#' functions return NA if the average distance or direction is not valid...
#' e.g., when averaging directions of 0 & 180 degrees, the result could
#' theoretically be 90 or 270 but is practically neither.
#' 
#' @param direction a vector of directions given in degrees (0 - 360) if
#' \code{deg}==TRUE or in radians if \code{deg}==FALSE
#' @param distance a vector of distances associated with each direction
#' @param deg a boolean object defining if \code{direction} is in degrees
#' (TRUE) or radians (FALSE)
#' @return \code{circular.averaging} returns the average direction while
#' \code{vector.averaging} returns a list with 2 elements distance & direction
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com} & Lorena Falconi
#' \email{lorefalconi@@gmail.com}
#' @examples
#' 
#' #EXAMPLE circular.averaging
#' circular.averaging(c(0,90,180,270)) #result is NA
#' circular.averaging(c(70,82,96,110,119,259))
#' 
#' #EXAMPLE vector.averaging
#' vector.averaging(c(10,20,70,78,108), distance=10)
#' vector.averaging(c(159,220,258,273,310),distance=runif(5))
#' 
#' @export 
circular.averaging = function(direction,deg=TRUE) {
	n=length(direction) #get the length of direction vector
	out = vector.averaging(direction=direction,distance=rep(1,n),deg=deg)
	return(out$direction)
}

#' @rdname circular.averaging
#' @export
vector.averaging = function(direction,distance,deg=TRUE) {
	if (deg) direction = direction*pi/180 #convert to radians
	n=length(direction) #get the length of direction vector
	if (any(is.na(direction))) { #ensure no NA data
		warning('NAs in data'); pos = which(is.na(direction)); direction = direction[-pos]; distance = distance[-pos] 
	} else {
		sinr <- sum(sin(direction))
		cosr <- sum(cos(direction))
		if (sqrt((sinr^2 + cosr^2))/n > .Machine$double.eps) {
			Ve = sum(distance*sin(direction))/n
			Vn = sum(distance*cos(direction))/n
			UV = sqrt(Ve^2 + Vn^2)
			AV1 = atan(Ve/Vn)
			AV2 = atan2(sinr, cosr)
			#perform some checks and correct when output in wrong quadrant
			AV = NULL
			if (abs(AV1-AV2) <= .Machine$double.eps) { #if both methods of determining directions are reporting the same values
				if (AV1>0) AV = AV1
				if (AV1<0) AV = 2 * pi + AV1		
			} else { #case when they are different values... add necessary values of pi
				if (AV1>0) AV = AV1 + pi
				if (AV1<0) AV = AV2
			}
			if (is.null(AV)) {
				return(list(distance=NA,direction=NA))
			} else {
				if (deg) AV = AV * 180 / pi #convert back to degrees
				return(list(distance=UV,direction=AV))
			}
		} else {
			return(list(distance=NA,direction=NA))
		}
	}
}
