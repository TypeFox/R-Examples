##' A method of the generic function rbind for objects of class trackdata
##' 
##' Different track data objects from one segment list are bound by combining
##' the \$data columns of the track data object by rows.  Track data objects
##' are created by \code{\link{get_trackdata}}.
##' 
##' All track data objects have to be track data of the same segment list.
##' Thus \$index and \$ftime values have to be identically for all track data
##' objects.  The number of columns of the track data objects must match. Thus
##' a track data object of more than one formant and single columned F0 track
##' data object can not be rbind()ed.
##' 
##' @aliases rbind.trackdata rbind
##' @param \dots track data objects
##' @return A track data object with the same \$index and \$ftime values of the
##' source track data objects and with \$data that includes all columns of
##' \$data of the source track data objects.
##' @author Jonathan Harrington
##' @seealso \code{\link{rbind}} \code{\link{cbind.trackdata}}
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
##'    vowlax.fund[1]
##'    
##'    #rms track data object for vowlax - first segment only 
##'    vowlax.rms[1]
##'    
##'    #now combine both track data objects
##'    fund.rms.lax = rbind(vowlax.fund[1:10,], vowlax.rms[1:10,])
##'  
##'    #the combined track data object
##'    #The first ten rows in \$data keep vowlax.fund data, the 11th to last row keeps vowlax.rms data 
##'    fund.rms.lax
##'    
##' 
##' 
##' @export rbind.trackdata
"rbind.trackdata" <- function (...) 
{
  mat <- NULL
  for (j in list(...)) {
    if (is.matrix(j$data)) 
      mat$data <- rbind(mat$data, j$data)
    else mat$data <- c(mat$data, j$data)
    mat$index <- rbind(mat$index, j$index)
    if (!is.null(j$ftime)) 
      mat$ftime <- rbind(mat$ftime, j$ftime)
  }
  diffinds <- mat$index[, 2] - mat$index[, 1] + 1
  right <- cumsum(diffinds)
  first.left <- diffinds - 1
  left <- right - first.left
  mat$index <- cbind(left, right)
  if (version$major >= 5) {
    oldClass(mat) <- "trackdata"
  }
  else {
    class(mat) <- "trackdata"
  }
  mat
}

