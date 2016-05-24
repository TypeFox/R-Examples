# $Id: elem.R 625 2005-06-09 14:20:30Z nj7w $

elem <- function(object=1, unit=c("KB","MB","bytes"), digits=0,
                 dimensions=FALSE)
{
  .Deprecated("ll", package="gdata")
  ll(pos=object, unit=unit, digits=digits, dimensions=dimensions)
}

