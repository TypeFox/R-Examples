clean.na <- function (x,margin,drop=FALSE)  {
  ## clean a matrix of all rows (margin=1) or columns (margin=2) with NA, 2000.02.23, C.Hoffmann
  ind <- apply(x,margin,function(xx) all(!is.na(xx)))
  if (margin == 1)   x[ind,,drop=drop]
  else x[,ind,drop=drop]
}

