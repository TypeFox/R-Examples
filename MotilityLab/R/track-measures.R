#' Track Measures
#'
#' Statistics that can be used to quantify tracks. All of these functions take a single
#' track as input and give a single number as output.
#'
#' @param x a single input track; a matrix whose first column is time and whose
#'  remaining columns are a spatial coordinate.
#' @param from index, or vector of indices, of the first row of the track. 
#' @param to index, or vector of indices, of last row of the track. 
#' @param xdiff row differences of x.
#'
#' @details
#' Some track measures consider only the first and last position (or steps) of a track, 
#' and are most useful in conjunction with \code{\link{aggregate.tracks}}; for instance,
#' \code{squareDisplacement} combined with \code{\link{aggregate.tracks}} gives a mean
#' square displacement plot, and \code{overallAngle} combined with 
#' \code{\link{aggregate.tracks}} gives a turning angle plot (see the examples for 
#' \code{\link{aggregate.tracks}}). To speed up computation of these measures on 
#' subtracks of the same track, the arguments \code{from}, \code{to} and 
#' possibly \code{xdiff} are exploited by \code{\link{aggregate.tracks}}.
#' 
#' \code{trackLength} sums up the distances between subsequent positsion; in other words,
#' it estimates the length of the underlying track by linear interpolation (usually
#' an underestimation). The estimation could be improved in some circumstances by using
#' \code{\link{interpolateTrack}}.
#'
#' \code{duration} returns the time elapsed between \code{x}'s first and last 
#' positions.
#'
#' \code{speed} simply divides \code{\link{trackLength}} by \code{\link{duration}}.
#' 
#' \code{displacement} returns the Euclidean distance between the track endpoints
#' and \code{squareDisplacement} returns the squared Euclidean distance.
#'
#' \code{displacementVector} returns the vector between the track endpoints.
#'
#' \code{maxDisplacement} computes the maximal Euclidean distance of any position 
#' on the track from the first position.
#'
#' \code{displacementRatio} divides the \code{displacement} by the \code{maxDisplacement}; 
#' \code{outreachRatio} divides the \code{maxDisplacement} by the \code{trackLength}
#' (Mokhtari et al, 2013). Both measures return 
#' values between 0 and 1, where 1 means a perfectly straight track.
#' If the track has \code{trackLength} 0, then \code{NaN} is returned.
#'
#' \code{straightness} divides the \code{displacement} by the \code{trackLength}.
#' This gives a number between 0 and 1, with 1 meaning a perfectly straight track. 
#' If the track has \code{trackLength} 0, then \code{NaN} is returned.
#'
#' \code{asphericity} is a different appraoch to measure straightness 
#' (Mokhtari et al, 2013): it computes the asphericity of the set of positions on the 
#' track _via_ the length of its principal components. Again this gives a number between 0 
#' and 1, with higher values indicating straighter tracks. 
#' Unlike \code{\link{straightness}}, however, asphericity ignores
#' back-and-forth motion of the object, so something that bounces between two positions 
#' will have low \code{straightness} but high \code{asphericity}. We define the 
#' asphericity of every track with two or fewer positions to be 1. For one-dimensional
#' tracks with one or more positions, \code{NA} is returned.
#'
#' \code{overallAngle} Computes the angle (in radians) between the first and the last 
#' segment of the given track. Angles are measured symmetrically, thus the return values
#' range from 0 to pi; for instance, both a 90 degrees left and right turns yield the
#' value pi/2. This function is useful to generate autocorrelation plots
#' (together with \code{\link{aggregate.tracks}}).
#' 
#' \code{meanTurningAngle} averages the \code{overallAngle} over all 
#' adjacent segments of a given track; a low \code{meanTurningAngle} indicates high
#' persistence of orientation, whereas for an uncorrelated random walk we expect 
#' 90 degrees. Note that angle measurements will yield \code{NA} values for tracks 
#' in which two subsequent positions are identical.
#'
#' \code{overallDot} computes the dot product between the first and the last 
#' segment of the given track. This function is useful to generate autocovariance plots
#' (together with \code{\link{aggregate.tracks}}).
#'
#' \code{hurstExponent} computes the corrected empirical Hurst exponent of the track.
#' This uses the function \code{\link[pracma]{hurstexp}} from the `pracma` package.
#' If the track has less than two positions, NA is returned.
#' \code{fractalDimension} estimates the fractal dimension of a track using the function
#' \code{\link[fractaldim]{fd.estim.boxcount}} from the 
#' `fractaldim` package. For self-affine processes in \eqn{n} dimensions, 
#' fractal dimension and Hurst exponent 
#' are related by the formula \eqn{H=n+1-D}. 
#' For non-Brownian motion, however, this relationship
#' need not hold. Intuitively, while the Hurst exponent takes a global approach to the 
#' track's properties, fractal dimension is a local approach to the track's properties
#' (Gneiting and Schlather, 2004).
#'
#' @name TrackMeasures
#'
#' @examples
#' ## show a turning angle plot with error bars for the T cell data.
#' with( (aggregate(BCells,overallDot,FUN="mean.se",na.rm=TRUE)),{
#'   plot( mean ~ i, xlab="time step", 
#'   	ylab="turning angle (rad)", type="l" )
#'   segments( i, lower, y1=upper )
#' } )
#'
#' @references
#' Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
#' and Marc Thilo Figge (2013), Automated Characterization and 
#' Parameter--Free Classification of Cell Tracks Based on Local Migration 
#' Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808
#'
#' Tillmann Gneiting and Martin Schlather (2004), Stochastic Models That Separate Fractal
#' Dimension and the Hurst Effect. \emph{SIAM Review} \bold{46}(2), 269--282. 
#' doi:10.1137/S0036144501394387
NULL

#' @rdname TrackMeasures
trackLength <- function(x) {
	if (nrow(x) > 2) { 
		dif <- apply(x[,-1,drop=FALSE], 2, diff)
		return(sum(sqrt(apply(dif^2, 1, sum))))
	} else if (nrow(x) == 2) {
		# this case is necessary because of dimension dropping by 'apply'
		return(sqrt(sum((x[2,-1] - x[1,-1])^2)))
	} else if (nrow(x) == 1) {
		return(0)
	} else {
		return(NA)
	}
}

#' @rdname TrackMeasures
duration <- function(x) {
	if( nrow(x) < 2 ){
		return(0)
	}
	dur <- x[nrow(x), 1] - x[1,1]
	dur <- unname(dur)
	return(dur)
}

#' @rdname TrackMeasures
speed <- function(x) {
  trackLength(x) / duration(x)
}

#' @rdname TrackMeasures
displacement <- function( x, from=1, to=nrow(x) ) {
	if( nrow(x) < 2 ){
		return(0)
	}
	sqrt(.rowSums((x[to, -1, drop=FALSE] - x[from, -1, drop=FALSE])^2,
		length(from),ncol(x)-1))
}

#' @rdname TrackMeasures
squareDisplacement <- function(x, from=1, to=nrow(x)) {
	if( nrow(x) < 2 ){
		return(0)
	}
	.rowSums((x[to, -1, drop=FALSE] - x[from, -1, drop=FALSE])^2,
		length(from),ncol(x)-1)
}

#' @rdname TrackMeasures
displacementVector <- function(x) {
  ret <- x[nrow(x),-1] - x[1,-1]
  rownames(ret) <- NULL  
  return(as.vector(ret))
}

#' @rdname TrackMeasures
maxDisplacement <- function(x) {
	limits <- c(1,nrow(x))
	sqrt(max(rowSums(sweep(x[seq(limits[1],limits[2]),-1,drop=FALSE],
		2,x[limits[1],-1])^2)))
}

#' @rdname TrackMeasures
displacementRatio <- function(x) {
  dmax <- maxDisplacement(x)
  if (dmax > 0) {
    return(displacement(x) / dmax)
  } else {
    return(NaN)
  }
}

#' @rdname TrackMeasures
outreachRatio <- function(x) {
  l <- trackLength(x) 
  if (l > 0) {
    return(maxDisplacement(x) / l)
  } else {
    return(NaN)
  }
}

#' @rdname TrackMeasures
straightness <- function(x) {
  l <- trackLength(x)
  if (l > 0) {
    return(displacement(x) / l)
  } else {
    return(1)
  }
}

#' @rdname TrackMeasures
overallAngle <- function(x, from=1, to=nrow(x), xdiff=diff(x)) {
	r <- rep(0, length(from))
	ft <- from<(to-1)
	if( sum(ft)>0 ){
		a <- xdiff[from[ft],-1,drop=FALSE]
		b <- xdiff[to[ft]-1,-1,drop=FALSE]
		a <- a/sqrt(.rowSums(a^2, nrow(a), ncol(a)))
		b <- b/sqrt(.rowSums(b^2, nrow(b), ncol(b)))
		rs <- .rowSums(a * b, nrow(a), ncol(a))
		rs[rs>1] <- 1
		rs[rs< -1] <- -1
		r[ft] <- acos(rs)
	}
	r
}

#' @rdname TrackMeasures
meanTurningAngle <- function(x) {
	if(nrow(x)<2){
		return(NaN)
	}
	mean(sapply(subtracks(x, 2), overallAngle),na.rm=TRUE)
}

#' @rdname TrackMeasures
overallDot <- function(x, from=1, to=nrow(x), xdiff=diff(x)) {
	r <- rep(NaN, length(from))
	ft <- from<to
	if( sum(ft) > 0 ){
		a <- xdiff[from[ft],-1,drop=FALSE]
		b <- xdiff[to[ft]-1,-1,drop=FALSE]
		r[ft] <- .rowSums(a * b, nrow(a), ncol(a))
	}
	r
}

#' @rdname TrackMeasures
asphericity <- function(x) {
	dim <- ncol(x) - 1
	limits <- c(1,nrow(x))
	if (limits[2]-limits[1]<2) {
		return(1)
	}
	if( dim == 1 ){
		return(NaN)
	}
	eigen.values <- eigen(cov(x[limits[1]:limits[2],-1]))$values
	rav <- mean(eigen.values)
	res <- sum((eigen.values - rav )^2 / dim / (dim - 1) / rav^2)
	return(res)
}

#' @rdname TrackMeasures
hurstExponent <- function(x) {
	if( !requireNamespace("pracma",quietly=TRUE) ){
		stop("This function requires the 'pracma' package.")
	}
  if (nrow(x) < 2) {
    return(NA)
  }
  return(pracma::hurstexp(diff(x[,-1,drop=FALSE]), display=FALSE)$Hal)
}

#' @rdname TrackMeasures
fractalDimension <- function(x){
  if( !requireNamespace("fractaldim",quietly=TRUE) ){
    stop("This function requires the 'fractaldim' package.")
  }
  return(fractaldim::fd.estim.boxcount(x[,-1])$fd)
}
