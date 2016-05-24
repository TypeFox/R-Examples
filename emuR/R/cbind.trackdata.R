##' A method of the generic function cbind for objects of class \'trackdata\'
##' 
##' Different track data objects from one segment list are bound by combining
##' the \$data columns of the track data object by columns.
##' 
##' All track data objects have to be track data of the same segment list.
##' Thus \$index and \$ftime values have to be identically for all track data
##' objects.  Track data objects are created by get_trackdata().  The number of
##' rows of the track data objects must match.
##' 
##' @aliases cbind.trackdata cbind
##' @param \dots track data objects
##' @return A track data object with the same \$index and \$ftime values of the
##' source track data objects and with \$data that includes all columns of
##' \$data of the source track data objects.
##' @author Jonathan Harrington
##' @seealso \code{\link{cbind}}, \code{\link{rbind.trackdata}}
##' \code{\link{trackdata}} \code{\link{get_trackdata}}
##' @keywords methods
##' @examples
##' 
##'    data(vowlax)
##'    
##'    #segment list vowlax - first segment only 
##'    vowlax[1,]
##'    
##'    #F0 track data object for vowlax - first segment only 
##'    vowlax.fund[1,]
##'    
##'    #rms track data object for vowlax - first segment only 
##'    vowlax.rms[1,]
##'       
##'    
##'    #now combine both track data objects
##'    fund.rms.lax = cbind(vowlax.fund, vowlax.rms)
##'    
##'    #the combined track data object - first segment only
##'    #The first column keeps vowlax.fund data, the second keeps vowlax.rms data 
##'    fund.rms.lax[1,]
##' 
##' 
##' @export cbind.trackdata
cbind.trackdata <- function (...) 
{
  mat <- NULL
  k <- 1
  for (j in list(...)) {
    if(k==1)
    {
      inds <- mat$index <- j$index
      mat$ftime <- j$ftime
    }
    else
    { 
      if( nrow(j$index) != nrow(inds) )
        stop("can't column bind trackdata from different segment lists")
      lvec = (j$index[,1]==inds[,1]) & (j$index[,2]==inds[,2])
      if(any(!lvec))
        stop("can't column bind trackdata from different segment lists")
    }
    k = k+1
  }
  
  
  
  for (j in list(...)) {
    mat$data <- cbind(mat$data, j$data)
  }
  
  if (version$major >= 5) {
    oldClass(mat) <- "trackdata"
  }
  else {
    class(mat) <- "trackdata"
  }
  mat
}
