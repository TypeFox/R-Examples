##' Smooth the data in a trackdata object.
##' 
##' Smooths each dataset in a trackdata object using a running mean smoother.
##' 
##' This function uses the \code{dapply} function to apply \code{smooth} to the
##' data for each segment.
##' 
##' @aliases dsmooth dsmooth.sub
##' @param dataset A trackdata object as returned from \code{track}.
##' @return The result of applying the \code{smooth} function to each column of
##' the data for each segment in the trackdata object.
##' @seealso smooth, dapply
##' @keywords misc
##' @export dsmooth
"dsmooth"<- function(dataset)
{
  ## dataset: a list, as returned by track
  ## separately smooths the segments corresponding to the data
  dapply(dataset, dsmooth.sub)
}



##' @export
"dsmooth.sub" <- function(data, ftime)
{
  if(is.matrix(data)){
    if(nrow(data)>5)
      data <- apply(data, 2, smooth)
  } else {
    if(length(data)>5)
      data <- smooth(data)
  }
  return( list(data=data, ftime=ftime) )
}
