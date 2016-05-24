#' Tracks Objects
#'
#' The function \code{tracks} is used to create tracks objects. \code{as.tracks} coerces
#' its argument to a tracks object, and \code{is.tracks} tests for tracks objects.
#' \code{c} can be used to combine (concatenate) tracks objects.
#'
#' @param x an object to be coerced or tested.
#'
#' @param ... for \code{tracks}, numeric matrices or objects that can be coerced to 
#'  numeric matrices. Each
#'  matrix contains the data of one track. The first column is the time, and the remaining
#'  columns define a spatial position. Every given matrix has to contain the same number
#'  of columns, and at least two columns are necessary.
#'
#'  For \code{c}, tracks objects to be combined.
#'
#'  For \code{as.tracks}, further arguments passed to methods (currently not used).
#'
#' @details Tracks objects are lists of matrices. Each matrix contains at least two 
#' columns; the first column is time, and the remaining columns are a spatial coordinate.
#' The following naming conventions are used (and enforced by \code{tracks}): The time
#' column has the name `t`, and spatial coordinate columns have names `x`,`y`,`z` if there
#' are three or less coordinates, and `x1`,...,`xk` if there are \eqn{k \ge 4} 
#' coordinates. All tracks in an object must have the same number of dimensions. The
#' positions in a track are expected to be sorted by time (and the constructor 
#' \code{tracks} enforces this).
#'
#' @examples
#' ## A single 1D track
#' x <- tracks( matrix(c(0, 8,  
#' 10, 9, 
#' 20, 7, 
#' 30, 7,
#' 40, 6, 
#' 50, 5), ncol=2, byrow=TRUE ) )
#'
#' ## Three 3D tracks
#' x2 <- tracks( rbind(
#'  c(0,5,0), c(1,5,3), c(2,1,3), c(3,5,6) ),
#'  rbind( c(0,1,1),c(1,1,4),c(2,5,4),c(3,5,1),c(4,-3,1) ),
#'  rbind( c(0,7,0),c(1,7,2),c(2,7,4),c(3,7,7) ) )
#'
#' @name tracks
NULL

#' Convert Tracks to Data Frame
#' 
#' Converts tracks from the list-of-matrices format, which is good
#' for efficient processing and therefore the default in this package, to a 
#' single dataframe which is convenient for plotting or saving the data.
#' 
#' @param x the \code{tracks} object to be coerced to a data frame.
#' 
#' @param row.names NULL or a character vector giving row names for the
#'  data frame.  Missing values are not allowed.
#' @param optional logical. Required for S3 consistency, but 
#' has no effect: column names are always assigned to the resulting
#'  data frame regardless of the setting of this option.
#' @param include.timepoint.column logical. If set to \code{TRUE}, then the resulting
#'  dataframe will contain a column that consecutively numbers the positions according
#'  to their time. Note that this information is anyway implicitly present in the time 
#'  information.
#' @param ... further arguments to be passed from or to other methods.
#'
#' @return A single data frame containing all individual tracks from the input with a 
#' prepended column named "id" containing each track's identifier in `x`.
#'
#' @examples
#' ## Display overall average position of the T cell data
#' colMeans( as.data.frame( TCells )[-c(1,2)] )
as.data.frame.tracks <- function(x, row.names = NULL, optional = FALSE, 
	include.timepoint.column=FALSE, ...) {
	ids <- rep(names(x), lapply(x,nrow))
	if( include.timepoint.column ){
		timepoint <- ave( ids, ids, FUN=seq_along )
		r <- data.frame(id=ids, timepoint=timepoint, do.call( 
			rbind.data.frame, x ) )
	} else {
		r <- data.frame(id=ids, do.call( 
			rbind.data.frame, x ) )
	}
	if( !is.null(row.names) && length(row.names)==nrow(r) ){
		rownames(r) <- row.names
	} else {
		rownames(r) <- NULL
	}
	return(r)
}

`[.tracks` <- function(x,y) as.tracks(as.list(x)[y])

#' @rdname tracks
as.tracks <- function(x, ...) 
  UseMethod("as.tracks")

#' @rdname tracks
as.tracks.list <- function(x, ...) 
  structure(x, class="tracks")


#' Convert from Tracks to List
#' 
#' Coerces a \code{tracks} object to a list.
#' 
#' @param x the \code{tracks} object to be coerced to a list.
#' @param ... further arguments to be passed from or to other methods.
as.list.tracks <- function(x, ...) 
  structure( x, class="list" )


#' @rdname tracks
is.tracks <- function(x)
  inherits(x, "tracks")

#' Filter Tracks 
#'
#' Extracts subtracks based on a given function.
#'
#' @param f a function that accepts a single track as its first argument and returns a 
#' logical value (or a value that can be coerced to a locical).
#' @param x a tracks object.
#' @param ... further arguments to be passed on to \code{f}.
#' @return A \code{tracks} object containing only those tracks from \code{x} for which
#' \code{f} evaluates to \code{TRUE}.
#' @examples
#' ## Remove short tracks from the T cells data
#' plot( filterTracks( function(t) nrow(t)>10, TCells ) )
filterTracks <- function(f,x,...){
	structure(x[as.logical(sapply(x,function(x) f(x,...) ))], class="tracks")
}

#' Sort Track Positions by Time
#' 
#' Sorts the positions in each track in a \emph{tracks} object by time.
#' 
#' @param x the \emph{tracks} object whose tracks are to be sorted by time.
#' 
#' @param decreasing logical.  Should the sort be increasing or decreasing? 
#'  Provided only for consistency with the generic sort method. The positions in
#'  each track should be sorted in increasing time order.
#' @param ... further arguments to be passed on to \code{order}.
#'
#' @details Sorts the positions of each track (represented as a data frame) in the 
#' \emph{tracks} object by time (given in the column t). 
#' 
#' @return A \emph{tracks object} that contains the tracks from the input object 
#' sorted by time is returned.
sort.tracks <- function(x, decreasing=FALSE, ...) {
	as.tracks(lapply(x, function(d) d[order(d[,"t"],decreasing=decreasing),,drop=FALSE] ))
}

#' @rdname tracks
c.tracks <- function(...) {
	args <- lapply(list(...), as.list)
	as.tracks(do.call(c, args))
}


#' Read Tracks from Text File
#' 
#' Reads cell tracks from a CSV or other text file. Data are expected to be organized as
#' follows.
#' One column contains a track identifier, which can be numeric or a string, and 
#' determines which points belong to the same track. 
#' Another column is expected to contain a time index or a time period (e.g. number of
#' seconds elapsed since the beginning of the track, or since the beginning of the 
#' experiment). Input of dates is not (yet) supported, as absolute time information is
#' frequently not available. 
#' One to three further columns contain the spatial coordinates 
#' (depending on whether the tracks are 1D, 2D or 3D). 
#' The names or indices of these columns in the CSV files are given using the 
#' corresponding parameters (see below). Names and indices can be mixed, e.g. you can
#' specify \code{id.column="Parent"} and \code{pos.columns=1:3}
#' 
#' @param file the name of the file which the data are to be read from, a 
#' readable text-mode connection or a complete URL
#' (see \code{\link[utils]{read.table}}).
#'
#' @param id.column index or name of the column that contains the track ID.
#'
#' @param time.column index or name of the column that contains elapsed time.
#'
#' @param pos.columns vector containing indices or names of the columns that contain
#' the spatial coordinates. If this vector has two entries and the second entry is NA,
#' e.g. \code{c('x',NA)} or \code{c(5,NA)} then all columns from the indicated column 
#' to the last column are used. This is useful when reading files where the exact number
#' of spatial dimensions is not known beforehand.
#'
#' @param scale.t a value by which to multiply each time point. Useful for changing units,
#' or for specifying the time between positions if this is not contained in the file
#' itself.
#'
#' @param scale.pos a value, or a vector of values, by which to multiply each spatial 
#' position. Useful for changing units.
#' 
#' @param header a logical value indicating whether the file contains the
#' names of the variables as its first line. See \code{\link[utils]{read.table}}.
#'
#' @param sep a character specifying how the colums of the data are separated. 
#' The default value \code{""} means columns are separated by tabs or other spaces.  
#' See \code{\link[utils]{read.table}}.
#'
#' @param track.sep.blankline logical. If set to \code{TRUE}, then tracks are expected
#' to be separated by one or more blank lines in the input file instead of being 
#' designated by a track ID column. In this case, numerical track IDs are automatically
#' generated.
#'
#' @param ... further arguments to be passed to \code{read.csv}, for instance 
#' \code{sep="\t"} can be useful for tab-separated files.
#' 
#' @return An object of class \emph{tracks} is returned, which is a list of 
#' matrices, each containing the positions of one track. The matrices 
#' have a column \eqn{t}, followed by one column for each of the input track's 
#' coordinates.
#' 
#' @details The input file's first four fields are interpreted as \eqn{id}, 
#' \eqn{pos}, \eqn{t} and \eqn{x}, respectively, and, if available, the fifth
#' as \eqn{y} and the sixth as \eqn{z}. The returned object has the class 
#' \emph{tracks}, which is a list of data frames representing the single 
#' tracks and having columns \eqn{t} and \eqn{x}, plus \eqn{y} and \eqn{z}, if
#' necessary. The tracks' ids are retained in their position in the list, while
#' the field \eqn{pos} will be unmaintained.
read.tracks.csv <- function(file, id.column=1, time.column=2, 
	pos.columns=c(3,4,5), scale.t=1, scale.pos=1, 
	header=TRUE, sep="", track.sep.blankline=FALSE, ...) {
	if( track.sep.blankline ){
		data.raw <- read.table(file, header=header, sep=sep,
			blank.lines.skip = FALSE, ...)
		is.blank <- apply( data.raw, 1, function(x) all( is.na(x) ) )
		ids <- cumsum( is.blank )
	} else {
		data.raw <- read.table(file, header=header, sep=sep, ...)
	}
	if( ncol(data.raw) < length(pos.columns)+1+as.integer(!track.sep.blankline) ){
		stop("CSV file does not contain enough columns! (Perhaps you need to specify 'sep')")
	}
	if( length(pos.columns) < 1 ){
		stop("At least one position column needs to be specified!")
	}
	if( length(pos.columns) == 2 && !is.finite(pos.columns[2]) ){
		cx <- match( pos.columns[1], colnames( data.raw ) )
		if( is.na(cx) && is.numeric(pos.columns[1]) ){
			cx <- pos.columns[1]
		}
		pos.columns <- seq( cx, length(data.raw) )
	}
	if( track.sep.blankline ){
		data.raw <- cbind( data.raw, ids )
		data.raw <- data.raw[!is.blank,]
		id.column <- ncol( data.raw )
	}
	cx <- as.character(c(id.column,time.column,pos.columns))
	cxc <- match( cx, colnames(data.raw) )
	cxi <- match( cx, seq_len(ncol(data.raw)) )

	cxi[is.na(cxi)] <- cxc[is.na(cxi)]

	if( any(is.na(cxi)) ){
		stop("Column(s) not found: ",
			paste(cx[is.na(cxi)],collapse=","))
	}
	r <- data.raw[,as.integer(cxi)]
	if( ncol(r) <= 5 ){
		colnames(r) <- c("id","t",c("x","y","z")[seq_along(pos.columns)])
	} else {
		colnames(r) <- c("id","t",paste0("x",seq_along(pos.columns)))
	}
	if( scale.t != 1 ){
		r[,"t"] <- scale.t*r[,"t"]
	}
	if( any( scale.pos != 1 ) ){
		r[,-c(1,2)] <- scale.t*r[,-c(1,2)]
	}
	sort.tracks(structure( 
		split.data.frame( as.matrix(r[,2:ncol(r)]), r[,1] ),
		class="tracks"))
}

#' Split Track into Multiple Tracks
#' 
#' @param x the input track (a data frame or a matrix).
#' @param id a string used to identify the resulting tracks; for instance,
#' if \code{id="X"}, then the resulting tracks are named X_1, X_2 and so forth.
#' Otherwise, they are simply labelled with integer numbers.
#' @param positions a vector of positive integers, given in ascending order.
#' @param min.length nonnegative integer. Resulting tracks that have fewer positions than
#'  the value of this parameter are dropped. 
splitTrack <- function( x, positions, id=NULL, min.length=2 ){
	freqs <- diff(c(0,positions,nrow(x)))
	segs <- rep( seq_len(length(positions)+1), freqs )
	segs[rep(freqs,freqs) < min.length] <- NA
	r <- structure( 
		split.data.frame( as.matrix(x), segs ),
		class="tracks")
	if(!is.null(id) && ( length(r) > 0 ) ){
		if( length(r) > 1 ){
			names( r ) <- paste0( id, "_", seq_along(r) )
		} else if( length(r) == 1 ){
			names( r ) <- id
		}
	}
	return(r)
}

#' Plot Tracks in 2D
#' 
#' Plots tracks contained in a "tracks" object into a twodimensional space
#' pallelel to the data's axes.
#' 
#' @param x the tracks to be plotted.
#' @param dims a vector giving the dimensions of the track data that shall be 
#' plotted, e.g. \code{c('x','y')} for the \eqn{x} and \eqn{y} dimension.
#' @param add boolean value indicating whether the tracks are to be added to the
#' current plot.
#' @param col a specification of the color(s) to be used. This can be a vector
#'  of size \code{length(x)}, where each entry specififes the color for the 
#'  corresponding track.
#' @param pch.start point symbol with which to label the first position of the track
#'  (see \code{\link[graphics]{points}}).
#' @param pch.end point symbol with which to label the last position of the track
#' @param cex point size for positions on the tracks.
#' @param ... additional parameters (e.g. xlab, ylab).
#' to be passed to \code{\link[graphics]{plot}}
#' (for \code{add=FALSE}) or \code{\link[graphics]{points}} (for \code{add=TRUE}), 
#' respectively.
#' @details One dimension of the data (by default \eqn{y}) is plotted against 
#' another (by default \eqn{x}). The dimesions can be chosen by means of the 
#' parameter \code{dims} and the axes can be labeled accordingly with the aid 
#' of \code{xlab} and \code{ylab}. The color can be set through \code{col}.
#' If the tracks should be added to an existing plot, \code{add} is to be set 
#' to \code{TRUE}.
#' @seealso \code{\link{plot3d}}
plot.tracks <- function(x, dims=c('x','y'), add=F,
		col=order(names(x)), pch.start=1, pch.end=NULL, 
		cex=.5, ... ) {
	args <- list(...)

	ids <- names(x)
	lcol <- rep_len(col, length(ids))
	pcol <- rep(lcol, lapply(x,nrow))
	dim1 <- unlist( lapply(x,'[', , dims[1]) )
	dim2 <- unlist( lapply(x,'[', , dims[2]) )

	if (add==T) {
		points( dim1, dim2, col=pcol, cex=cex, ... )
	} else {
		plot( dim1, dim2, col=pcol, cex=cex, ...)
	}
  
	if( !is.null(pch.start) ){
		starting.points <- lapply( x, function(p) p[1,] )
		px <- sapply(starting.points,'[[',dims[1])
		py <- sapply(starting.points,'[[',dims[2])
		points( px, py, col=col, pch=pch.start, cex=2 )
	}

	if( !is.null(pch.end) ){
		end.points <- lapply( x, function(p) p[nrow(p),] )
		px <- sapply(end.points,'[[',dims[1])
		py <- sapply(end.points,'[[',dims[2])
		points( px, py, col=col, pch=pch.end, cex=2 )
	}
  
	lseg <- function(f,i) {
		function(d) f( d[,dims[i]], -1 ) 
	}

	x0 <- .ulapply(x,lseg(head,1))
	y0 <- .ulapply(x,lseg(head,2))

	x1 <- .ulapply(x,lseg(tail,1))
	y1 <- .ulapply(x,lseg(tail,2))

	lcol <- rep(lcol, lapply(x, function(d) {
		nrow(d) - 1
	}))
	segments(x0, y0, x1, y1, col=lcol)
}

#' Create Track Object from Single Track
#' 
#' Makes a \code{tracks} object containing the given track.
#' 
#' @param x the input track.
#' 
#' @return A list of class \code{tracks} containing only the input track \code{x}, which
#' is assigned the name "1". 
wrapTrack <- function(x) {
  return(as.tracks(list("1" = x)))
}

#' Staggered Version of a Function
#' 
#' Returns the "staggered" version of a track measure. That is, instead of
#' computing the measure on the whole track, the measure is averaged over
#' all subtracks (of any length) of the track.
#' 
#' @param measure a track measure (see \link{TrackMeasures}).
#' @param ... further parameters passed on to \code{\link{applyStaggered}}.
#' 
#' @details This is a wrapper mainly designed to provide a convenient interface
#' for track-based staggered computations with \code{lapply}, see example.
#' 
#' @return Returns a function that computes
#' the given measure in a staggered fashion on that track.
#' 
#' @references 
#' Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
#' and Marc Thilo Figge (2013), Automated Characterization and 
#' Parameter--Free Classification of Cell Tracks Based on Local Migration 
#' Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808
#'
#' @examples
#' hist( sapply( TCells, staggered( displacement ) ) )
#'
staggered <- function(measure, ...){
	return( function(x) applyStaggered(x, measure, ...) )
}


#' Compute a Measure on a Track in a Staggered Fashion
#' 
#' Computes a measure on all subtracks of a track and return them either
#' as a matrix or return their mean.
#' 
#' @param x the track for which the measure is to be computed.
#' @param measure the measure that is to be computed.
#' @param matrix a logical indicating whether the whole matrix of values for 
#' the measure for each of the input track's subtracks is to be returned. 
#' Otherwise only the mean is returned.
#' @param min.segments the number of segments that each regarded subtrack 
#' should at least consist of. Typically, this value would be set to the 
#' minimum number of segments that a (sub)track must have in order for the 
#' measure to be decently computed. For example, at least two segments are needed
#' to compute the \code{\link{overallAngle}}.
#' 
#' @details The measure is computed for each of the input track's subtracks of 
#' length at least \code{min.segments}, and the resulting values are either 
#' returned in a matrix (if \code{matrix} is set), or their mean is returned. 
#' The computed matrix is symmetric since the direction along which a 
#' track is traversed is assumed not to matter. The values at 
#' \code{[i, i + j]}, where j is a nonnegative integer with 
#' \eqn{j < }\code{min.segments}, (with the default value \code{min.segments=1}
#' this is exactly the main diagonal) are set to \code{NA}.
#'  
#' @return If \code{matrix} is set, a matrix with the values of the measure for
#' all the input track's subtracks is returned. The value of this matrix at 
#' position \code{[i, j]} corresponds to the subtrack that starts with the input track's
#' \eqn{j}th point and ends at its \eqn{i}th. Thus, with increasing column number,
#' the regarded subtrack's starting point is advanced on the original track, 
#' and the values for increasingly long subtracks starting from this point can 
#' be found columnwise below the main diagonal, respectively.
#' If `matrix` is not set, the mean over the values of the measure for all
#' subtracks of at least `min.segments` segments is retruned.
#' 
#' @examples 
#' ## Compute the staggered matrix for overallAngle applied to all long enough 
#' ## subtracks of the first T cell track
#' applyStaggered(TCells[[1]], overallAngle, matrix=TRUE, min.segments = 2)
applyStaggered <- function(x, measure, matrix=FALSE, min.segments=1) {
	if (matrix) {
		mat <- matrix(nrow=nrow(x), ncol=nrow(x), 0)
		diag(mat) <- NA
	} else {
		stag.meas <- c()
	}
	if (!matrix) {
		segMeans <- .computeSegmentwiseMeans(x, measure, min.segments)
		n.before <- length(segMeans)
		segMeans <- segMeans[!is.nan(segMeans)]
		n <- length(segMeans)
		if (n.before!=n) {
		  warning("NaNs have been removed.")
		}
		weights <- seq(n, 1)
		weights <- weights[!is.nan(segMeans)]
		ret <- weighted.mean(segMeans, weights)
		return(ret)
	} else {
		for (i in (min.segments):(nrow(x) - 2)) {
			subtracks <- subtracks(x, i)
			val <- sapply(subtracks, measure)
			diag(mat[-1:-i,]) <- val
		}
		mat[nrow(x), 1] <- sapply(subtracks(x, (nrow(x)-1)), measure)
		mat <- mat + t(mat)
		return(mat)
	}
}

#' Decompose Track(s) into Subtracks
#' 
#' Creates a \emph{tracks} object consisting of all subtracks of `x`
#' with `i` segments (i.e., `i`+1 positions).
#' 
#' @param x a single track or a \code{tracks} object.
#' @param i subtrack length. A single integer, lists are not supported.
#' @param overlap the number of segments in which each subtrack shall overlap 
#' with the previous and next subtrack. The default \code{i - 1} returns all
#' subtracks. Can be negative, which means that space will be left between
#' subtracks. 
#' 
#' @details The output is always a single \emph{tracks} object, which is 
#' convenient for many common analyses. If subtracks are to be considered separately
#' for each track, use the function \code{\link{staggered}} together with
#' \code{lapply}. Subtrack extraction always starts at the first position of the
#' input track.
#' 
#' @return A \emph{tracks} object is returned which contains all the subtracks 
#' of any track in the input \emph{tracks} object that consist of exactly `i`
#' segments and overlap adjacent subtracks in `overlap` segments.
subtracks <- function(x, i, overlap=i-1 ) {
	if( !is.tracks(x) ){
		x <- wrapTrack( x )
	}
	do.call(c, lapply(x, 
  		function(t) .computeSubtracksOfISegments(t, i, overlap )))
}

#' Get Track Prefixes 
#' 
#' Creates a \code{tracks} object consisting of all prefixes (i.e., subtracks
#' starting with the first position of a track) of `x`
#' with `i` segments (i.e., `i`+1 positions).
#' 
#' @param x a single track or a \code{tracks} object.
#' @param i subtrack length. A single integer, lists are not supported.
#'
#' @details This function behaves exactly like \code{\link{subtracks}} except
#' that only subtracks starting from the first position are considered.
#'
prefixes <- function(x,i) {
	if( !is.tracks(x) ){
		x <- wrapTrack( x )
	}
	lapply( x,
		function(t) {
			if( nrow(t) <= i ){ return( NULL ) }
			t[1:(i+1), ,drop=FALSE]
		}
	)
}

#' Normalize a Measure to Track Duration
#' 
#' Returns a measure that divides the input measure by the duration of its 
#' input track.
#' 
#' @param x a track measure (see \link{TrackMeasures}).
#' 
#' @return A function that computes the input measure for a given track
#' and returns the result divided by the track's duration.
#' 
#' @examples 
#' ## normalizeToDuration(displacement) can be used as an indicator 
#' ## for the motion's efficiency
#' sapply(TCells, normalizeToDuration(displacement))
normalizeToDuration <- function(x) {
  function(track) {
    x(track) / duration(track)
  }
}



#' Bivariate Scatterplot of Track Measures
#' 
#' Plots the values of two measures applied on the given tracks against each 
#' other.
#' 
#' @param measure.x the measure to be shown on the X axis (see \link{TrackMeasures}).
#' @param measure.y the measure to be shown on the Y axis.
#' @param x the input \code{tracks} object.
#' @param add a logical indicating whether the tracks are to be added to an 
#' existing plot via \code{\link[graphics]{points}}.
#' @param xlab label of the x-axis. By default the name of the input function 
#' \code{measure.x}.
#' @param ylab label of the y-axis. By default the name of the input function 
#' \code{measure.y}.
#' @param ... additional parameters to be passed to \code{\link[graphics]{plot}} 
#' (in case \code{add=FALSE}) or \code{\link[graphics]{points}} (\code{add=TRUE}).
# For ellipse documentation:
#' @inheritParams hotellingsTest
#' 
#' @details Plots the value of \code{measurey} applied to \code{x} against the
#' value of \code{measurey} applied to \code{y}. This is useful for "FACS-like"
#' motility analysis, where clusters of cell tracks are identified based on their
#' motility parameters (Moreau et al, 2012; Textor et al, 2014).
#'
#' @references
#' Moreau HD, Lemaitre F, Terriac E, Azar G, Piel M, Lennon-Dumenil AM,
#' Bousso P (2012), Dynamic In Situ Cytometry Uncovers 
#' T Cell Receptor Signaling during Immunological Synapses and Kinapses In Vivo.
#' \emph{Immunity} \bold{37}(2), 351--363. doi:10.1016/j.immuni.2012.05.014
#'
#' Johannes Textor, Sarah E. Henrickson, Judith N. Mandl, Ulrich H. von Andrian,
#' J\"urgen Westermann, Rob J. de Boer and Joost B. Beltman (2014),
#' Random Migration and Signal Integration Promote Rapid and Robust T Cell Recruitment.
#' \emph{PLoS Computational Biology} \bold{10}(8), e1003752.
#' doi:10.1371/journal.pcbi.1003752
#'
#' @examples
#' ## Compare speed and straightness of 3 example population tracks.
#' ## To make the comparison fair, analyze subtracks of fixed length.
#' plotTrackMeasures( subtracks(TCells,4,0), speed, straightness, ellipse.col="black" )
#' plotTrackMeasures( subtracks(BCells,4,0), speed, straightness,
#'   col=2, ellipse.col=2, pch=2, add=TRUE )
#' plotTrackMeasures( subtracks(Neutrophils,4,0), speed, straightness, 
#'   col=3, ellipse.col=3, pch=3, add=TRUE )
plotTrackMeasures <- function(x, measure.x, measure.y, add=FALSE, 
                              xlab=deparse(substitute(measure.x)), 
                              ylab=deparse(substitute(measure.y)), 
                              ellipse.col="red", ellipse.border="black",
                              conf.level=0.95, ...) {
	if( !is.tracks(x) ){
		x <- as.tracks( x )
	}
	mx <- sapply(x,measure.x)
	my <- sapply(x,measure.y)
	if(add) {
		points(mx, my, ...)
	} else {
		plot(sapply(x,measure.x), sapply(x,measure.y), xlab=xlab, ylab=ylab, ...)
	}
	if( !is.na(ellipse.col) ){
		n <- length(mx)
		p <- 1
		df.1 <- p
		df.2 <- n-p
		d <- cbind( mx, my )
		S <- cov( d )
		t <- sqrt(((n-1)*df.1/(n*df.2))*qf(1-conf.level,df.1,df.2,lower.tail=F))
		polygon( ellipse::ellipse( S, centre=colMeans(d),
			t=t ),
			border=ellipse.border,
			col=.setColAlpha(ellipse.col,128) )
	}
}

#' Cluster Tracks
#' 
#' Perform a hierarchical clustering of a set of tracks according to a given vector
#' of track measures.
#' 
#' @param tracks the tracks that are to be clustered.
#' @param measures a function, or a vector of functions (see \link{TrackMeasures}). 
#' Each function is expected to 
#' return a single number given a single track.
#' @param scale logical indicating whether the measures values shall be scaled
#' using the function \code{\link[base]{scale}} before the clustering. 
#' @param ... additional parameters to be passed to \code{\link[stats]{hclust}}.
#'
#' @return An object of class *hclust*, see \code{\link[stats]{hclust}}.
#' 
#' @details The measures are applied to each of the tracks in the given
#' \emph{tracks} object. According to the resulting values, the tracks are 
#' clustered using a hierarchical clustering (see \code{\link[stats]{hclust}}).
#' If \code{scale} is \code{TRUE}, the measure values are scaled to mean value 
#' \eqn{0} and standard deviation \eqn{1} (per measure) before the clustering.
#'
#' @examples 
#' ## Cluster tracks according to the mean of their Hust exponents along X and Y
#'
#' cells <- c(TCells,Neutrophils)
#' real.celltype <- rep(c("T","N"),c(length(TCells),length(Neutrophils)))
#' ## Prefix each track ID with its cell class to evaluate the clustering visually
#' names(cells) <- paste0(real.celltype,seq_along(cells))
#' clust <- clusterTracks( cells, hurstExponent )
#' plot( clust )
#' ## How many cells are "correctly" clustered?
#' sum( real.celltype == c("T","N")[cutree(clust,2)] )
#' 
clusterTracks <- function(tracks, measures, scale=TRUE, ...) { 
	values <- matrix(nrow=length(tracks))
	if (is.function(measures)) {
		measures <- c(measures)
	}
	values <- do.call(cbind, lapply(measures, function(m) sapply(tracks, m)))
	if (scale) {
		values <- scale(values)
	}
	hclust(dist(values), ...)
}

#' Compute Summary Statistics of Subtracks
#' 
#' Computes a given measure on subtracks of a given track set, applies a summary 
#' statistic for each subtrack length, and returns the results in a convenient form.
#' This important workhorse function facilitates many common motility analyses 
#' such as mean square displacement, turning angle, and autocorrelation plots.
#' 
#' @param x the tracks object whose subtracks are to be considered. 
#' If a single track is given, it will be coerced to a tracks object 
#' using \code{\link{wrapTrack}} (but note that this requires an explicit call
#' \code{aggregate.tracks}).
#' @param measure the measure that is to be computed on the subtracks.
#' @param by a string that indicates how grouping is performed. Currently, two
#' kinds of grouping are supported: 
#' \itemize{
#'  \item{"subtracks"}{ Apply \code{measure} to all subtracks according to 
#' the parameters \code{subtrack.length} and \code{max.overlap}.}
#'  \item{"prefixes"}{ Apply \code{measure} to all prefixes (i.e., subtracks starting
#'  from a track's initial position) according to the parameter \code{subtrack.length}.}
#' }
#' 
#' @param FUN a summary statistic to be computed on the measures of subtracks with the
#' same length. Can be a function or a string.
#' If a string is given, it is first matched to the following builtin values:
#' \itemize{
#'  \item{"mean.sd"}{ Outputs the mean and \eqn{mean - sd} as lower and 
#'  \eqn{mean + sd} as upper bound}
#'  \item{"mean.se"}{ Outputs the mean and \eqn{mean - se} as lower and 
#'  \eqn{mean + se} as upper bound}
#'  \item{"mean.ci.95"}{ Outputs the mean and upper and lower bound of a 
#'  parametric 95 percent confidence intervall}. 
#'  \item{"mean.ci.99"}{ Outputs the mean and upper and lower bound of a 
#'  parametric 95 percent confidence intervall}. 
#'  \item{"iqr"}{ Outputs the interquartile range, that is, the median, and the 
#'  25-percent-quartile as a lower and and the 75-percent-quartile as an 
#'  upper bound}
#' }
#' If the string is not equal to any of these, it is passed on to 
#' \code{\link{match.fun}}.
#'
#' @param subtrack.length an integer or a vector of integers defining which subtrack
#' lengths are considered. In particular, \code{subtrack.length=2} 
#' corresponds to a "step-based analysis" (Beltman et al, 2009).
#'
#' @param max.overlap an integer controlling what to do with overlapping subtracks.
#' A maximum overlap of \code{max(subtrack.length)} will imply
#' that all subtracks are considered. For a maximum overlap of 0, only non-overlapping
#' subtracks are considered. A negative overlap can be used to ensure that only subtracks
#' a certain distance apart are considered. In general, for non-Brownian motion there will
#' be correlations between subsequent steps, such that a negative overlap may be necessary
#' to get a proper error estimate.
#' @param na.rm logical. If \code{TRUE}, then \code{NA}'s and \code{NaN}'s 
#' are removed prior to computing the summary statistic.
#' @param filter.subtracks a function that can be supplied to exclude certain subtracks
#' from an analysis. For instance, one may wish to compute angles only between steps of
#' a certain minimum length (see Examples).
#' @param ... further arguments passed to or used by methods.
#'
#' @details For every number of segments \eqn{i} in the set defined by 
#' \code{subtrack.length}, all subtracks of any track in the input 
#' \code{tracks} object that consist of exactly \eqn{i} segments are 
#' considered. The input \code{measure} is applied to the subtracks individually, 
#' and the \code{statistic} is applied to the resulting values. 
#' 
#' @return A data frame with one row for every \eqn{i} 
#' specified by \code{subtrack.length}. The first column contains the values 
#' of \eqn{i} and the remaining columns contain the values of the summary statistic
#' of the measure values of tracks having exactly \eqn{i} segments.
#' @examples 
#' ## A mean square displacement plot with error bars.
#' dat <- aggregate(TCells, squareDisplacement, FUN="mean.se")
#' with( dat ,{
#'   plot( mean ~ i, xlab="time step", 
#'   	ylab="mean square displacement", type="l" )
#'   segments( i, lower, y1=upper )
#' } )
#'
#' ## Compute a turning angle plot for the B cell data, taking only steps of at least
#' ## 1 micrometer length into account
#' check <- function(x) all( sapply( list(head(x,2),tail(x,2)), trackLength ) >= 1.0 )
#' plot( aggregate( BCells, overallAngle, subtrack.length=1:10, 
#'   filter.subtracks=check )[,2], type='l' )
#'
#' ## Compare 3 different variants of a mean displacement plot
#' # 1. average over all subtracks
#' plot( aggregate( TCells, displacement ), type='l' )
#' # 2. average over all non-overlapping subtracks
#' lines( aggregate( TCells, displacement, max.overlap=0 ), col=2 )
#' # 3. average over all subtracks starting at 1st position
#' lines( aggregate( TCells, displacement, by="prefixes" ), col=3 )
#'
#' @references
#' Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009),
#' Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
#' 789--798. doi:10.1038/nri2638
aggregate.tracks <- function( x, measure, by="subtracks", FUN=mean, 
    subtrack.length=seq(1, (maxTrackLength(x)-1)),
    max.overlap=max(subtrack.length), na.rm=FALSE, 
    filter.subtracks=NULL, ... ){
    if( class( x ) != "tracks" ){
    	if( class( x ) %in% c("data.frame","matrix" ) ){
    		x <- wrapTrack( x )
    	} else {
    		stop("Cannot coerce argument 'tracks' to tracks object" )
    	}
    }
    if (max(subtrack.length) > (maxTrackLength(x)-1)) {
		warning("No track is long enough!")
  	}
  	the.statistic <- NULL
  	if (is.character(FUN)) {
		if (FUN == "mean.ci.95") {
		  the.statistic <- function(x) {
		  	mx <- mean(x,na.rm=na.rm)
			ci <- tryCatch( t.test(x)$conf.int, error=function(e) rep(NA,2) )
			return(c(lower=ci[1], mean=mx, upper=ci[2]))
		  }
		} else if (FUN == "mean.ci.99") {
		  the.statistic <- function(x) {
			mx <- mean(x,na.rm=na.rm)
			ci <- tryCatch( t.test(x, conf.level=.99)$conf.int, 
				error=function(e) rep(NA,2) )
			return(c(lower=ci[1], mean=mean(x), upper=ci[2]))
		  }
		} else if (FUN == "mean.se") {
		  the.statistic <- function(x) {
		  	mx <- mean(x,na.rm=na.rm)
		  	sem <- sd(x,na.rm=na.rm)/sqrt(sum(!is.na(x))) 
		  		# note that this also works is na.omit is FALSE
			return(c(mean=mx, lower = mx - sem, upper = mx + sem))
		  }
		} else if (FUN == "mean.sd") {
		  the.statistic <- function(x) {
		  	mx <- mean(x,na.rm=na.rm)
		  	sem <- sd(x,na.rm=na.rm)
			return(c(mean=mx, lower = mx - sem, upper = mx + sem))
		  }
		} else if (FUN == "iqr") {
		  the.statistic <- function(x) {
			lx <- quantile(x,probs=c(.25,.5,.75),na.rm=na.rm)
			names(lx) <- c("lower","median","upper")
			return(lx)
		  }
		} else {
			the.statistic.raw <- match.fun( FUN )
		}
	} else {
		if ( !is.function( FUN) ) {
		   stop("FUN must be a function or a string")
		} else {
			the.statistic.raw <- FUN
		}
	}
	if( is.null( the.statistic ) ){
		if( na.rm ){
			the.statistic <- function(x){
				the.statistic.raw(x[!is.na(x) & !is.nan(x)])
			}
		} else {
			the.statistic <- the.statistic.raw
		}
	}
	if( is.function( filter.subtracks ) ){
		isValidSubtrack <- function(t) !is.null(t) && filter.subtracks(t)
	} else {
		isValidSubtrack <- Negate(is.null)
	}
	if( by == "subtracks" ){
		# optimized version: avoids making copies of subtracks (relevant especially
		# for measures like displacement, which actually only consider the track
		# endpoints). Also pre-computes the diff of the track if needed.
		if( ( length( intersect(c("x","limits"),names(formals(measure))) ) == 2 )
			&& !is.function( filter.subtracks ) ){
			if( "xdiff" %in% names(formals(measure)) ){
				ft <- function(t,i) apply( .subtrackIndices(i,min(max.overlap,i-1),nrow(t)), 
					1, measure, x=t, xdiff=diff(t) )
			} else {
				ft <- function(t,i) apply( .subtrackIndices(i,min(max.overlap,i-1),nrow(t)), 
					1, measure, x=t )
			}
			measure.values <- list()
			k <- 1
			xi <- x
			for( i in sort(subtrack.length) ){
				xi <- xi[sapply(xi,nrow)>i]
				measure.values[[k]] <- .ulapply( xi, ft, i=i )
				k <- k+1
			}
		} else if( ( length( intersect(c("x","from","to"),names(formals(measure))) ) == 3 )
			&& !is.function( filter.subtracks ) ){
			measure.values <- list()
			k <- 1
			xi <- x
			for( i in sort(subtrack.length) ){
				xi <- xi[sapply(xi,nrow)>i]
				xi1 <- do.call( rbind, xi )
				si <- .subtrackIndices(i,min(max.overlap,i-1),sapply(xi,nrow))
				measure.values[[k]] <- measure(xi1,si[,1],si[,2])
				k <- k+1
			}
		} else {
			# unoptimized version: makes copies of all subtracks.
			the.subtracks <- list()
			measure.values <- list()
			k <- 1
			for (i in subtrack.length) {
				the.subtracks[[k]] <- subtracks(x, i, min(max.overlap,i-1) )
				the.subtracks[[k]] <- Filter( isValidSubtrack, the.subtracks[[k]] )
				measure.values[[k]] <- sapply( the.subtracks[[k]], measure )
				k <- k+1
			}
		}
	} else if (by == "prefixes" ) {
		measure.values <- lapply( subtrack.length, 
				function(i) sapply( Filter( isValidSubtrack, prefixes( x, i )  ), 
					measure ) 
			)
	} else {
		stop( "grouping method ", by, " not supported" )
	}
	# this is necessary because of automatic insertion of NULL elements when e.g.
	# measure.values[[5]] is assigned to
	measure.values <- Filter( Negate(is.null), measure.values )
	value <- sapply(measure.values, the.statistic)
	ret <- rbind(i=subtrack.length, value)
	return(data.frame(t(ret)))
}

#' Length of Longest Track
#' 
#' Determines the maximum number of positions over the tracks in \code{x}.
#' 
#' @param x the \code{tracks} object the tracks in which are to be considered.
#' 
#' @return The maximum number of rows of a track in \code{x}
maxTrackLength <- function(x) {
  return(max(sapply(x, nrow)))
}


#' Bounding Box of a Tracks Object
#' 
#' Computes the minimum and maximum coordinates per dimension (including time) 
#' for all positions in a given list of tracks.
#' 
#' @param x the input \code{tracks} object.
#'  
#' @return Returns a matrix with two rows and \eqn{d+1} columns, where \eqn{d} is 
#' the number of spatial dimensions of the tracks. The first row contains the minimum 
#' and the second row the maximum value of any track in the dimension given by 
#' the column.
#'
#' @examples
#' ## Use bounding box to set up plot window
#' bb <- boundingBox(c(TCells,BCells,Neutrophils))
#' plot( Neutrophils, xlim=bb[,"x"], ylim=bb[,"y"], col=1 )
#' plot( BCells, col=2, add=TRUE )
#' plot( TCells, col=3, add=TRUE )
boundingBox <- function(x) {
	if( !is.tracks(x) ){
		x <- as.tracks(x) 
	}
	bounding.boxes <- lapply(x, .boundingBoxTrack) 
	empty <- mat.or.vec(nrow(bounding.boxes[[1]]), ncol(bounding.boxes[[1]]))
	colnames(empty) <- colnames(bounding.boxes[[1]])
	rownames(empty) <- rownames(bounding.boxes[[1]])
	empty[1,] <- Inf
	empty[2,] <- -Inf
	return(Reduce(.minMaxOfTwoMatrices, bounding.boxes, empty)) 
}

#' Test Unbiasedness of Motion
#' 
#' Test the null hypothesis that a given set of tracks originates from an uncorrelated
#' and unbiased type of motion (e.g., a random walk without drift). This is done by
#' testing whether the mean step vector is equal to the null vector. 
#' 
#' @param tracks the tracks whose biasedness is to be determined.
#' @param dim vector with the names of the track's 
#' dimensions that are to be considered. By default c("x", "y").
#' @param plot logical indicating whether the scatter of the step's directions, 
#' origin of ordinates (green circle) and the mean of the data points (green 
#' cross) are to be plotted. (In one dimension also the bounds of the 
#' condfidence interval are given.) Plot works only in one or two dimensions.
#' @param add whether to add the plot to the current plot (\code{TRUE}) or create a 
#" new one (\code{FALSE}).
#' @param step.spacing How many positions are to be left out between 
#' the steps that are considered for the test. For persistent motion, subsequent
#' steps will be correlated, which leads to too low p-values because Hotelling's 
#' test assumes that the input data is independent. To avoid this, the resulting
#' p-value should either be corrected for this dependence (e.g. by adjusting 
#' the degrees of freedom accordingly), or `step.spacing` should be set to a value
#' high enough to ensure that the considered steps are approximately independent.
#' @param ellipse.col color with which to draw the confidence ellipse of the mean (for
#'  1D, this corresponds to the confidence interval of the mean). 
#'  Use \code{NA} to omit the confidence ellipse.
#' @param ellipse.border color of the confidence ellipse border. Use \code{NA} to omit
#'  the border.
#' @param conf.level the desired confidence level for the confidence ellipse.
#' @param ... further arguments passed on to \code{plot}.
#' 
#' @return A list with class \code{htest}.
#' 
#' @details Computes the displacement vectors of all segments in the tracks 
#' given in \code{tracks}, and performs Hotelling's T-square Test on that vector.
#' 
#' @examples 
#' ## Test H_0: T-cells migrate by uncorrelated random walk on x and y coordinates,
#' ## and report the p-value.
#' hotellingsTest( TCells )$p.value
#'
#' @references
#' Johannes Textor, Antonio Peixoto, Sarah E. Henrickson, Mathieu
#'  Sinn, Ulrich H. von Andrian and Juergen Westermann (2011),
#'	Defining the Quantitative Limits of Intravital Two-Photon Lymphocyte Tracking.
#' \emph{PNAS} \bold{108}(30):12401--12406. doi:10.1073/pnas.1102288108
#'
hotellingsTest <- function(tracks, dim=c("x", "y"), 
	step.spacing=0, plot=FALSE, add=FALSE, ellipse.col="blue", 
	ellipse.border="black", conf.level=0.95, ... ) {
	if( plot && !(length(dim) %in% c(1,2) ) ){
		stop("Plotting is only supported for 1D and 2D")
	}
	Tx <- projectDimensions(tracks, dim)
	sx <- sapply( subtracks(Tx, 1, overlap=-step.spacing), displacementVector )
	if( length(dim) > 1 ){
		sx <- t(sx)
	} else {
		sx <- matrix(sx,ncol=1)
	}
	n <- dim(sx)[1]
	p <- dim(sx)[2]
	df.1 <- p
	df.2 <- n-p
	S <- cov(sx)
	mX <- colMeans(sx)
	A <- solve(S)
	T2 <- n * t(mX) %*% A %*% mX
	p.value <- 1-pf(T2*df.2/df.1/(n-1), df.1, df.2)
	if (plot && (length(dim)==2)) {
		if (plot) {
			if( add ){
				points(sx,...)
			} else {
				plot(sx,...)
			}
			abline(h=0)
			abline(v=0)
			if( !is.na(ellipse.col) ){
				t <- sqrt(((n-1)*df.1/(n*df.2))*qf(1-conf.level,df.1,df.2,lower.tail=F))
				polygon(ellipse::ellipse(S,centre=mX,t=t),
					col=.setColAlpha(ellipse.col,128),border=ellipse.border)
				points(mX[1],mX[2],col=ellipse.col,pch=20)
			}
		}
	} else if (plot && (length(dim)==1)) {
		sd.sx <- sd(sx)
		plot(0, 0, col=3, xlab = dim[1], ylim=c(-1,1), ...)
		stripchart(sx, add=TRUE, at=0, pch=1)
		abline(v=0)
		t <- sqrt(((n-1)*df.1/(n*df.2))*qf(1-conf.level,df.1,df.2,lower.tail=F))
		polygon(ellipse::ellipse(rbind(c(S,0),c(0,.1)),centre=c(mX,0),t=t),
			col=.setColAlpha(ellipse.col,128),border=NA)
		points(mX[1],0,col=ellipse.col,pch=20)
	}
	structure(list(statistic=c(T2=T2),
		p.value=p.value,
		method="Hotelling's one sample T2-test",
		parameter=c(df1=df.1,df2=df.2),
		alternative="two.sided",
		null.value=c(location=paste0("c(", paste0(rep(0,ncol(sx)), collapse = ","), ")"))
		),
		class="htest")
}

#' Select Tracks by Measure Values
#' 
#' Given a tracks object, extract a subset based on upper and lower bounds of a certain
#' measure. For instance, extract all tracks with a certain minimum length.
#'
#' @param x the input tracks.
#' @param measure measure on which the selection is based (see \link{TrackMeasures}).
#' @param lower specifies the lower bound (inclusive) of the allowable measure.
#' @param upper specifies the upper bound (inclusive) of the allowable measure.
#' @examples
#' ## Slower half of T cells
#' slow.tcells <- selectTracks( TCells, speed, -Inf, median( sapply(TCells,speed) ) )
selectTracks <- function(x,measure,lower,upper){
	if(is.tracks(x)){
		mdat <- sapply(x,measure)
		return(x[mdat<=upper & mdat>=lower])
	}else{
		obj <-deparse(substitute(tracks))
		warning(sprintf("Obj:%s is not a tracks object",obj))
	}
}

#' Plot Tracks in 3D
#'
#' Takes an input tracks object and plots them in 3D using the 
#' \link[scatterplot3d]{scatterplot3d} function.
#'
#' @param x the tracks which will be plotted in 3d
#' @param ... further arguments to be passed on to 
#' \link[scatterplot3d]{scatterplot3d}
#'
#' @examples
#' if( require("scatterplot3d",quietly=TRUE) ){
#'   plot3d( TCells )
#' }
plot3d <- function(x,...){
	if( !requireNamespace("scatterplot3d",quietly=TRUE) ){
		stop("This function requires the package 'scatterplot3d'.")
	}
	tracks_df <- as.data.frame.tracks(
		lapply( x, function(t) rbind(t,rep(NA,ncol(t)) ) ) )
	s3d <- scatterplot3d::scatterplot3d(tracks_df[,-c(1,2)],
		type="n",xlab="X Position",ylab="Y Position",zlab="Z Position",...)
	colvec <- rainbow(length(names(x)))[tracks_df[,1]]
	pts <- s3d$xyz.convert( tracks_df[,-c(1,2)] )
	segments( head(pts$x,-1), head(pts$y,-1), tail(pts$x,-1), 
		tail(pts$y,-1), col=colvec )
}

#' Compute Time Step of Tracks
#'
#' Applies a summary statistics on the time intervals between pairs of consecutive 
#' positions in a track dataset. 
#' 
#' @param x the input tracks.
#' @param FUN the summary statistic to be applied.
#' @param na.rm logical, indicates whether to remove missing values before applying FUN.
#'
#' @details Most track quantification depends on the assumption that track positions are
#' recorded at constant time intervals. If this is not the case, then most of the statistics
#' in this package (except for some very simple ones like \code{\link{duration}}) will not work.
#' In reality, at least small fluctuations of the time steps can be expected. This
#' function provides a means for quality control with respect to the tracking time.
#'
#' @return Summary statistic of the time intervals between two consecutive positions in a track 
#' dataset. 
#'
#' @examples
#' ## Show tracking time fluctuations for the T cell data
#' d <- timeStep( TCells )
#' plot( sapply( subtracks( TCells, 1 ) , duration ) - d, ylim=c(-d,d) ) 
timeStep <- function( x, FUN=median, na.rm=FALSE ){
	if( class(FUN) == "character" ){
		FUN=match.fun( FUN )
	}
	x <- .ulapply( x, function(x) diff(x[,1]) )
	if( na.rm ){
		x <- x[!is.na(x)]
	}
	return( FUN( x ) )
}

#' Interpolate Track Positions
#'
#' Approximates the track positions at given time points using linear interpolation
#' (via the \code{\link[stats]{approx}} function).
#'
#' @param x the input track (a matrix or data frame).
#' @param t the times at which to approximate track positions. These must lie 
#' within the interval spanned by the track timepoints.
#' @param how specifies how to perform the interpolation. Possible values are 
#' \code{"linear"} (which uses \code{\link[stats]{approx}} with default values) and
#' \code{"spline"} (which uses \code{\link[stats]{spline}} with default values). 
#'
#' @examples
#' ## Compare interpolated and non-interpolated versions of a track
#' bb <- boundingBox( TCells[2] )
#' plot( TCells[2] )
#' t2i <- interpolateTrack(TCells[[2]], seq(bb[1,"t"],bb[2,"t"],length.out=100),"spline")
#' plot( tracks( t2i ), add=TRUE, col=2 )
interpolateTrack <- function( x, t, how="linear" ){
	f=approx
	if( how=="spline" ){
		f=spline
	}
	pos.interpolated <- apply(x[,-1],2,function(p) f(x[,1], p, xout=t)$y) 
	r <- cbind( t, pos.interpolated )
	colnames(r) <- colnames(x)
	return(r)
}

#' @rdname tracks
tracks <- function( ... ){
	tlist <- lapply( list(...), data.matrix )
	if( !all( sapply( tlist , is.numeric ) ) ){
		stop("Track time/position information must be numeric!")
	}
	tdims <- sapply( tlist, ncol )
	if( any( tdims < 2 ) ){
		stop("At least time and one spatial dimension are needed!")
	}
	if( length(tlist) == 0 ){
		return( as.tracks(tlist) ) 
	}
	if( length(tlist) > 0 &&
		diff( range( tdims ) ) > 0 ){
		stop("All tracks must have the same number of dimensions!")
	}
	r <- tlist[[1]]
	if( ncol(r) <= 5 ){
		cls <- c("t",c("x","y","z")[seq_len(ncol(r)-1)])
		for( i in seq_along(tlist) ){
			colnames(tlist[[i]]) <- cls
		}
	} else {
		cls <- c("t",paste0("x",seq_len(ncol(r)-1)))
		for( i in seq_along(tlist) ){
			colnames(tlist[[i]]) <- cls
		}
	}
	if( is.null(names(tlist)) ){
		names(tlist) <- seq_along(tlist)
	}
	sort(as.tracks(tlist))
}

#' Subsample Track by Constant Factor
#'
#' Make tracks more coarse-grained by keeping only every \emph{k}th position.
#'
#' @param x an input track or tracks object.
#' @param k a positive integer. Every \eqn{k}th position of each
#' input track is kept.
#' @seealso \code{interpolateTrack}, which can be used for more flexible track 
#' coarse-graining. 
#' @examples
#' ## Compare original and subsampled versions of the T cell tracks
#' plot( TCells, col=1 )
#' plot( subsample( TCells, 3 ), col=2, add=TRUE, pch.start=NULL )
subsample <- function( x, k=2 ){
	if( is.tracks(x) ){
		as.tracks( lapply( x, function(t) subsample(t,k) ) )
	} else {
		x[seq(1,nrow(x),by=k),]
	}
}
