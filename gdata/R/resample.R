## The S/R 'sample' function behaves differently if it is passed a
## sampling vector of length 1 than if it is passed a
## vector of length greater than 1.  For the 1-element
## case it samples from the list 1:x, instead of from the contents
## of x.  This function remove the special case: it always samples from
## the provided argument, no matter the length.
resample <- function(x, size, replace = FALSE, prob = NULL)
  {
    if(length(x)<1)
      if(!missing(size) && size>0)
        stop("Requested sample of size ", size, " from list of length 0")
      else
        x[FALSE]
    else if(length(x)==1)
      {
        if(missing(size) || size==1)
          x
        else if(size>=1 && replace==TRUE)
          rep(x, size)
        else if(size < 1)
          x[FALSE]
        else
          stop("Cannot cannot take a sample larger than the population",
               " when 'replace = FALSE'")
      }
    else
      sample(x, size, replace, prob)
  }
