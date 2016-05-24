##' Get the track times from EMU trackdata objects
##' 
##' The function obtains the times at which track values occur.
##' 
##' Every \$data value in a trackdata object is associated with a time at which
##' it occurs in the utterance. This function returns those times.
##' 
##' @param trackdata An EMU trackdata object, or a matrix of track values
##' obtained at a single time point using dcut()
##' @author Jonathan Harrington
##' @seealso \code{\link{start.trackdata}} \code{\link{end.trackdata}}
##' \code{\link{start.emusegs}} \code{\link{end.emusegs}}
##' @keywords datagen
##' @examples
##' 
##' # track time values for a trackdata object
##' times <- tracktimes(vowlax.fdat)
##' # track time values for a matrix of trackdata values
##' # at  the temporal midpoint
##' tracktimes(dcut(vowlax.fdat[1:3,], 0.5, prop=TRUE))
##' 
##' @export tracktimes
"tracktimes" <- function(trackdata)
{
  if(is.trackdata(trackdata))
    # return the times at which the frames
    # of trackdata occur as a numerical vector
    times <- as.numeric(dimnames(trackdata$data)[[1]])
  else if(is.vector(trackdata))
    times <- as.numeric(names(trackdata))
  else if(is.matrix(trackdata))
    times <- as.numeric(dimnames(trackdata)[[1]])
  else times <- NULL
  times
}
