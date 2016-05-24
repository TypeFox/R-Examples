##' Extract a subset of data from a trackdata object
##' 
##' A function that cuts up trackdata either at a proportional time or
##' proportionally between two times. It is a subsidiary function of dplot()
##' 
##' 
##' @aliases dextract dextract.sub
##' @param dataset A trackdata object
##' @param start A single valued numeric vector corresponding to a proportional
##' time between zero (the onset of the trackdata) and one (the offset of the
##' trackdata).
##' @param end As start, but optional
##' @return If both start and end are specified, a trackdata object is
##' returned, otherwise a vector if the original trackdata is one-dimensional
##' and the end argument is not used, or a matrix if the original trackdata has
##' more than one dimension and the end argument is not used
##' @author Jonathan Harrington
##' @seealso \code{dcut}
##' @keywords datagen
##' @examples
##' 
##' data(demo.vowels.f0)
##' data(demo.vowels.fm)
##' 
##' form = demo.vowels.fm
##' # get the formants at the midpoint: f50 is a matrix
##' # same as dcut(form, .5, prop=TRUE)
##' f50 = dextract(form, 0.5)
##' # get the formants between the 25% and 75% time points
##' # fcut is a trackdata object
##' # same as dcut(form, .25, .75, prop=TRUE)
##' fcut = dextract(form, 0.25, 0.75)
##' # get  F0 at the midpoint. fzero50 is a vector
##' # same as dcut(fzero, .5, prop=TRUE)
##' fzero = demo.vowels.f0
##' fzero50 = dextract(fzero, 0.5)
##' 
##' 
##' @export dextract
dextract <- function(dataset, start, end) {
  if((start < 0) | (start > 1)) {
    stop("proportional duration must be between 0 and 1")
  }
  if(!missing(end)) {
    if((end < 0) | (end > 1)) {
      stop("proportional duration must be between 0 and 1")
    }
    if(start > end) {
      stop("proportional start must be less than proportional end")
    }
  }
  
  if(missing(end)) {
    leftin <- dataset$index[, 1]
    rightin <- dataset$index[, 2]
    scalein <- round((rightin - leftin) * start)
    outin <- leftin + scalein
    if(is.matrix(dataset$data))
      return(dataset$data[outin,  ])
    else return(dataset$data[outin])
  } else {
    dapply(dataset, dextract.sub, start, end)
  }
}


##' @export
"dextract.sub" <- function (data, ftime, start, end) 
{
# helper function for use via dapply, returns a new
# trackdata element cut at start/end proportions
  len <- nrow(data)
  start <- floor(start * (len - 1) + 1)
  end <- ceiling(end * (len - 1) + 1)
  
  newdata <- data[start:end, ]
  times <- seq(ftime[1], ftime[2], length = len)
  newftime <- times[c(start, end)]
  return(list(data = newdata, ftime = newftime))
}
