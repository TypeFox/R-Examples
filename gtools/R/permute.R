# $Id: permute.R 625 2005-06-09 14:20:30Z nj7w $

permute <- function(x) sample( x, size=length(x), replace=FALSE )
