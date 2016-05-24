
##' A method of the generic function dim for objects of class 'trackdata'
##' 
##' The function returns the dimension attributes of a track data object.
##' 
##' The function returns the dimension attributes of a track data object as the
##' number of segments x number of tracks.  c(nrow(x$index), ncol(x$data))
##' 
##' @aliases dim.trackdata dim
##' @param x a track data object
##' @author Jonathan Harrington
##' @keywords methods
##' @examples
##' 
##'    #isol.fdat is the formant track of the segment list isol
##' 
##'    #write out the dimension of the track data object 
##'    dim(isol.fdat)
##' 
##'    #because there are 13 segments
##'    isol.fdat$ftime
##' 
##'    #and there are 4 rows for each segment (see here for the first segment)
##'    isol.fdat$data[1,]
##' 
##' @export
dim.trackdata <- function(x)
{
  # function returns the dimension attributes of
  # a trackdata object as the number of segments x number of tracks
  c(nrow(x$index), ncol(x$data))
}




##' Dimnames of trackdata object
##' 
##' returns dimension names of trackdata objects
##' 
##' 
##' @param x trackdata object
##' @keywords methods
##' @export
"dimnames.trackdata" <- function(x)
{
  trackdata = x
  dimnames(trackdata$data)
}
