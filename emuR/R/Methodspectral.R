##' Expand spectral
##' 
##' see function
##' 
##' 
##' @aliases [.spectral
##' @keywords internal
##' @export
"[.spectral" <- function (dat, i, j, drop) 
{
  
  if(!is.trackdata(dat))
  {
    if(is.matrix(dat))
    {
      if(missing(j))
        j <- freqtoint(dat, trackfreq(dat))
      else
      {
        if(is.logical(j))
          j <- trackfreq(dat)[j]
        j <- freqtoint(dat, j)
        
      }
      o <- NextMethod("[")
      class(o) <- c(class(o), "spectral")
      
      attr(o, "fs") <- attr(dat, "fs")[j]
      return(o)
      
    }
    else{
      
      if(missing(i))
        i <- freqtoint(dat, trackfreq(dat))
      else
      {
        if(is.logical(i))
          i <- trackfreq(dat)[i]
        i <- freqtoint(dat, i)
      }
      o <- NextMethod("[")
      class(o) <- c(class(o), "spectral")
      attr(o, "fs") <- attr(dat, "fs")[i]
      return(o)
    }
  }
}
