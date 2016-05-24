# midpoints, ninenum, seventeennum
# http://code.google.com/p/cowares-excel-hello/source/browse/trunk/ellipseplot/
#
# Copyright (C) 2013 Tomizono
# Fortitudinous, Free, Fair, http://cowares.nobody.jp
#

# return list of divided range of c(min.index, max.index)
midpoint <- function(x) {
  n <- (x[1] + x[2]) / 2
  list(c(x[1], floor(n)), c(ceiling(n), x[2]))
}

# dig midpoints of array
midpoints.calc <- function(x, n) {
  nx <- length(x)
  pl <- list(c(1, nx))
  for(i in 1:n) pl <- unlist(lapply(pl, midpoint), recursive=F)
  pm <- matrix(c(1, simplify2array(pl), nx), ncol=2, byrow=T)
  apply(pm, 1, function(a) sum(x[a]) / 2)
}

# midpoints
midpoints <- function(x, n=1, na.rm=TRUE) {
  xna <- is.na(x)
  if(na.rm) x <- x[!xna]
  if((!na.rm && any(xna)) || (length(x) == 0)) 
    return(rep.int(NA, 2^n + 1))

  midpoints.calc(sort(x), n)
}

# derived
ninenum <- function(x, na.rm=TRUE) midpoints(x, 3, na.rm)
seventeennum <- function(x, na.rm=TRUE) midpoints(x, 4, na.rm)

