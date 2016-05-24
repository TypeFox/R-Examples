##' Build trackdata objects from the output of by()
##' 
##' buildtrack() converts a list that is the output of by.trackdata() into a
##' trackdata object if the list components are matrices whose rows are
##' successive values in time.
##' 
##' The default of by.trackdata() is to return a list. If each element of the
##' list consists of a matrix whose rows are values occurring at the times
##' given by the row dimension names of the matrix, then buildtrack() can be
##' used to convert the list into a trackdata object. If the times are not
##' given in the row dimension names, then these can be supplied as an
##' additional argument to buildtrack()
##' 
##' @param mylist a list that ist output from by()
##' @param ftime ftime
##' @param trackname name of track data object
##' @author Jonathan Harrington
##' @seealso \code{\link{by}}
##' @keywords manip
##' @examples
##' 
##'    #vowlax.fdat is a track data objects of formant of the vowlax segment list
##'    #calculate the difference between adjacent formant values
##'    p = by(vowlax.fdat[1,2],INDICES=NULL, diff)
##'    
##'    p
##'    
##'    
##'    #now build a track data object out of these values
##'    m = buildtrack(p)
##'    
##'    m
##' 
##' @export buildtrack
"buildtrack" <- function(mylist, ftime=NULL, trackname="")
{
  # convert a list that is usually output from by() into
  # a trackdata object
  
  # examples
  # p = by(vowlax.fdat[1:3,2], diff)
  # m = trackfromlist(p)
  # p = by(vowlax.fdat[1:3,], apply, 2, diff)
  # o = trackfromlist(p)
  # r = by(vowlax.fdat[1:3,2], smooth)
  # should break because no ftimes
  # trackfromlist(r)
  # e = trackfromlist(r, ftime=vowlax.fdat[1:3,]$ftime)
  # m = by(vowlax.fdat[1:3,], apply, 2, smooth)
  # e = trackfromlist(m, ftime=vowlax.fdat[1:3,]$ftime)
  
  res <- NULL
  for(j in 1:length(mylist)){
    if(!is.matrix(mylist[[j]]))
      mylist[[j]] <- cbind(mylist[[j]])
    N <- nrow(mylist[[j]])
    if(is.null(ftime))
    {
      if(is.null(dimnames(mylist[[j]])[[1]]))
        stop("can't find any ftime values")
      times <- as.numeric(dimnames(mylist[[j]])[[1]])
      times <- c(times[1], times[N])
    }
    else 
      times <- ftime[j,]
    res$data <- rbind(res$data, cbind(mylist[[j]]))
    res$index <- c(res$index, N)
    res$ftime <- rbind(res$ftime, times)
  }
  # build indices
  n <- res$index
  right <- cumsum(n)
  left <- c(1, right+1)
  left <- left[-length(left)]
  res$index <- cbind(left, right)
  as.trackdata(res$data, res$index, res$ftime, trackname)
}
