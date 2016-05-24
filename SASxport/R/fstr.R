fstr <- function(name, length, digits)
  {
    invalid <- function(x) is.null(x) | ( length(x)<1 ) | ( nchar(x) < 1 ) | x==0
    inner <- function(i)
      {
        if( invalid(name[i]) )
          return("")
        if( invalid( length[i] ) )
          return(name[i])
        if( invalid(digits[i]) )
          return( paste(name[i], length[i], '.', sep='' ) )
        else
          return( paste(name[i], length[i], '.', digits[i], sep='' ) )
      }
    if(length(name)>0)
      sapply( 1:length(name), inner)
    else
      character(0)
  }
