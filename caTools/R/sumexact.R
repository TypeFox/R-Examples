#===========================================================================#
# caTools - R library                                                       #
# Copyright (C) 2005 Jarek Tuszynski                                        #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

sumexact = function(..., na.rm = FALSE)
{
  x = c(...,  recursive=TRUE)
  if (na.rm) x = x[!is.na(x)]
  else if (any(is.na(x))) return(NA)
  n = length(x)
  .C("sum_exact", as.double(x), y = as.double(0), as.integer(n),
      NAOK=TRUE, DUP=TRUE, PACKAGE="caTools")$y
}

#==============================================================================

cumsumexact = function(x)
{
  n = length(x)
  .C("cumsum_exact", as.double(x), y = double(n), as.integer(n),
      NAOK=TRUE, DUP=TRUE, PACKAGE="caTools")$y
}


