

### makemovie ###


#' @export
#' @name makemovie
#' @aliases LSD.makemovie
#' @title Interpolate rows of a matrix to extend the number of cols
#' @description Interpolate rows of a matrix to yield a smooth transitions.
#' @param input a matrix with numerical entries.
#' @param timepoints a integer vector containing the timepoints corresponding to the columns of 'input'.
#' @param timestep a non-negative integer specifying the number of timesteps between the existing timepoints (defaults to \code{1}, if not specified).
#' @param motionline a integer vector giving the timepoints of the resulting matrix (derived from timepoints and timesteps by default).
#' @author Achim Tresch, Bjoern Schwalb
#' @seealso \code{\link{clusterplot}}, \code{\link{align}}, \code{\link{demotour}}
#' @examples len = 10
#' x = sin(seq(0,2*pi,length=len*2))
#' fun = function(){n=sample(1:len,1);return(x[n:(n+len-1)])}
#' input = t(replicate(7,fun(),simplify=TRUE))
#' input = input + rnorm(length(input))/2
#' par(mfrow=c(1,2))
#' plotmatrix(input,main="original",cols=1:7,type="o")
#' mov = makemovie(input,timestep=0.2)
#' plotmatrix(mov,main="interpolated",cols=1:7,type="o")
#' @keywords matrix


makemovie = function(input,timepoints=NULL,timestep=1,motionline=NULL)
{
	if (!is.matrix(input)) stop("'input' must be a matrix !")
	if (sum(is.na(input)) != 0) stop("'input' must not contain NAs !")
	if (is.null(timepoints)) timepoints = 1:ncol(input)
	sorted = order(timepoints,decreasing=FALSE)
	input = input[,sorted]
	timepoints = timepoints[sorted]
	colnames(input) = timepoints
	if (is.null(motionline)){
		nrpics = ceiling(diff(range(timepoints))/timestep)+1
		motionline = seq(timepoints[1],timepoints[length(timepoints)],length=nrpics)		
	}
	timeline = function(x){splinefun(timepoints,x)(motionline)}
	movie = t(apply(input,1,timeline))
	colnames(movie) = round(motionline,digits=2)
	return(movie)
}


### aliases ###


LSD.makemovie = makemovie



