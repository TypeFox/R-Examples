##' Extract a subset of data from a trackdata object
##' 
##' Extract a subset of data from a trackdata object according to the label in
##' the parallel label vector.
##' 
##' 
##' @param dataset A trackdata object returned from \code{track}.
##' @param labs A vector of labels parallel to \code{trackdata$index}, i.e. one
##' for each segment in the trackdata.
##' @param labtype A vector of labels for which data is to be extracted.
##' @return A trackdata object which is a subset of \code{trackdata} containing
##' only the data for those labels in \code{labtype}.  The result has the same
##' components as the input \code{trackdata}:
##' 
##' \item{data}{ A vector or matrix of numerical data. } \item{index}{ A two
##' column matrix giving the start and end indeces into the data vector for
##' each segment. } \item{ftime}{ A two column matrix giving the start and end
##' times for each segment. }
##' @seealso track, dextract, get.time.element, frames.time
##' @keywords internal
##' @export dextract.lab
"dextract.lab"<- function(dataset, labs, labtype = unique(labs))
{
  # extract data values from a dataset ($data, $index, $ftime)
  # according to labtype (e.g. "i:", c("i:", "u:").
  # labs is parallel to dataset$index; labtype are
  # the label types for which the values in dataset are
  # to be extracted
  mat <- NULL
  lvals <- dataset$index[, 2] - dataset$index[, 1] + 1
  newlabs <- rep(labs, lvals)
  temp <- muclass(newlabs, labtype)
  if(is.matrix(dataset$data))
    vals <- dataset$data[temp,  ]
  else 
    vals <- dataset$data[temp]
  
  temp.lab <- muclass(labs, labtype)
  
  if(!is.null(dataset$ftime))
    ftimes <- dataset$ftime[temp.lab,  ]
  
  finds <- dataset$index[temp.lab,  ]	
  ## readjust the indeces
  diffinds <- finds[, 2] - finds[, 1] + 1
  right <- cumsum(diffinds)
  first.left <- diffinds - 1
  left <- right - first.left
  finds <- cbind(left, right)
  mat$data <- vals
  mat$index <- finds
  if(!is.null(dataset$ftime))
    mat$ftime <- ftimes
  if( version$major >= 5 ) {
    oldClass(mat) <- "trackdata"
  } else {
    class(mat) <- "trackdata"
  }
  mat
}
