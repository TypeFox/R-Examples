# $Id: invalid.R 625 2005-06-09 14:20:30Z nj7w $

invalid <- function(x)
  {
    if( missing(x) || is.null(x) || length(x)==0 )
      return(TRUE)
    if(is.list(x))
      return(all(sapply(x,invalid)))
    else if(is.vector(x))
      return(all(is.na(x)))
    else
      return(FALSE)
  }
