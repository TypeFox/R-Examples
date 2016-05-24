
#' Normalize Tracks
#' 
#' Translates each track in a given set of tracks such that the 
#' first position is the origin. 
#' 
#' @param x the input \code{tracks} object.
#'
#' @examples
#' ## normalization of Neutrophil data reveals upward motion
#' plot( normalizeTracks( Neutrophils ) )
normalizeTracks <- function(x){
	as.tracks(lapply(x, .normalizeTrack))
}


#' Extract Spatial Dimensions
#' 
#' Projects tracks onto the given spatial dimensions.
#' 
#' @param x the input tracks object.
#' @param dims a character vector (for column names) or an integer vector (for column
#'  indices) giving the dimensions to extract from each track.
#'  The time dimension (i.e., the first column of all tracks) is always included.
#' 
#' @return A tracks object is returned that contains only those dimensions 
#' of the input \code{tracks} that are given in \code{dims}.
#'
#' @examples
#' ## Compare 2D and 3D speeds
#' speed.2D <- mean( sapply( subtracks( projectDimensions( TCells, c("x","z") ), 2 ), speed ) )
#' speed.3D <- mean( sapply( TCells, speed ) )
#' 
projectDimensions <- function(x, dims=c("x","y")) {
	if( class(x) != "tracks" ){
		x <- as.tracks(x)
	}
	if( length(x) == 0 ){
		return(x)
	}
	if( class(dims) == "character" ){
		.dims <- match(dims,colnames(x[[1]]))
		if( any(is.na(.dims)) ){
			stop("dimensions not found: ",dims[is.na(.dims)])
		}
		as.tracks(lapply(x, function(t) {
			t[, c(1,.dims)]
		}))
	} else if( class(dims) == "numeric" ){
		as.tracks(lapply(x, function(t) {
			t[, c(1,dims)]
		}))
	} else {
		stop("'dims' must be a character or integer vector!")
	}
}

#' Process Tracks Containing Gaps
#'
#' Many common motility analyses, such as mean square displacement plots, assume that
#' object positions are recorded at constant time intervals. For some application domains,
#' such as intravital imaging, this may not always be the case. This function can be 
#' used to pre-process data imaged at nonconstant intervals, provided the deviations are
#' not too extreme.
#'
#' @param x the input tracks object.
#' @param how string specifying what do with tracks that contain gaps. Possible
#'   values are:
#' \itemize{
#'  \item{"drop":}{ the simplest option -- discard all tracks that contain gaps.}
#'  \item{"split":}{ split tracks around the gaps, e.g. a track for which the step
#'  between the 3rd and 4th positions is too long or too short is split into one
#'  track corresponding to positions 1 to 3 and another track corresponding to
#'  position 3 onwards.}
#'  \item{"interpolate":}{ approximate the track positions using linear
#'  interpolation (see \code{\link{interpolateTrack}}). The result is a tracks
#'  object with constant step durations.
#'  }
#' }
#' @param tol nonnegative number specifying by which fraction each step may deviate 
#'  from the average step duration without being considered a gap. For instance, if
#'  the average step duration (see \code{\link{timeStep}}) is 100 seconds and \code{tol}
#'  is 0.05 (the default), then step durations between 95 and 105 seconds (both inclusive)
#'  are not considered gaps. This option is ignored for \code{how="interpolate"}.
#' @param split.min.length nonnegative integer. For \code{how="split"}, this
#' discards all resulting tracks shorter than
#' this many positions. 
#'
#' @examples
#' ## The Neutrophil data are imaged at rather nonconstant intervals
#' print( length( Neutrophils ) )
#' print( length( repairGaps( Neutrophils, tol=0.01 ) ) )
repairGaps <- function( x, how="split", tol=0.05, split.min.length=2 ){
	deltaT <- timeStep( x, na.rm=TRUE )
	if( how=="drop" ){
		gi <- which( sapply( x, function(t) length( .gaps( t, tol=tol, deltaT=deltaT ) ) >0 ) )
		if( length(gi)>0 ){
			return( as.tracks( x[-gi] ) )
		} else {
			return(x)
		}
	} else if( how=="split" ){
		ids <- names(x)
		as.tracks( unlist( lapply( seq_along(x), function(i){
			splitTrack( x[[i]], .gaps( x[[i]], tol, deltaT ), id=ids[i], 
				min.length=split.min.length )
		} ), recursive=FALSE ) )
	} else if( how=="interpolate" ){
		as.tracks( lapply( x, function(t)
			interpolateTrack( t, seq(t[1,1],t[nrow(t),1],by=deltaT) ) ) )
	} else {
		stop("Invalid value for parameter \"how\"")
	}
}

