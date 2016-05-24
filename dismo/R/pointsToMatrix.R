# Author: Robert J. Hijmans & Jacob van Etten
# October 2009
# version 1
# license GPL3


.pointsToMatrix <- function(p, checkLonLat=TRUE) {
	if (inherits(p, 'SpatialPoints')) {
		test <- !is.projected(p)
		if (! isTRUE (test) ) {
			if (is.na(test)) {
				warning('Coordinate reference system of SpatialPoints object is not set. Assuming it is degrees (longitude/latitude)!')  			
			} else {
				stop('Points are projected. They should be in degrees (longitude/latitude)')  
			}
			# or rather transform them ....?
		}
		p <- coordinates(p)
	}
	if (is.data.frame(p)) {
		p <- as.matrix(p)
	}
	if (is.vector(p)){
		if (length(p) != 2) {
			stop('Wrong length for a vector, should be 2')
		} else {
			p <- matrix(p, ncol=2) 
		}
	} else if (is.matrix(p)) {
		if (length(p[1,]) != 2) {
			stop( 'A points matrix should have 2 columns')
		}
		cn <- colnames(p)
		if (length(cn) == 2) {
			if (toupper(cn[1]) == 'Y' | toupper(cn[2]) == 'X')  {
				warning('Suspect column names (x and y reversed?)')
			}
			if (toupper(substr(cn[1],1,3) == 'LAT' | toupper(substr(cn[2],1,3)) == 'LON'))  {
				warning('Suspect column names (longitude and latitude reversed?)')
			}
		}		
	} else {
		stop('points should be vectors of length 2, matrices with 2 columns, or a SpatialPoints* object')
	}

	if (checkLonLat) {
		if (min(p[,1], na.rm=TRUE) < -720) { stop('longitude < -720') }
		if (max(p[,1], na.rm=TRUE) > 720) {  stop('longitude > 720')  }
		if (min(p[,1], na.rm=TRUE) < -360) { warning('longitude < -360') }
		if (max(p[,1], na.rm=TRUE) > 360) {  warning('longitude > 360')  }
		if (min(p[,2], na.rm=TRUE) < -90) {  stop('latitude < -90')  }
		if (max(p[,2], na.rm=TRUE) > 90) {  stop('latitude > 90')  }
	}
	
	if (! is.numeric(p) ) { p[] <- as.numeric(p) }
	
	return(p)
}

