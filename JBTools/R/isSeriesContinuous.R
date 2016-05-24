isSeriesContinuous = function(
    ##title<< Test for continuous (non NA interrupted) series
    x ##<< numeric vector: series to test
  )
  ##description<<
  ## isSeriesContinuous test whether a vector contains one non interrupted sequence
  ## of values (i.e. no NaN in between). 
  ##details<<
  ## The function returns TRUE when the vector contains one (and only one)
  ## non NA interupted sequence of values. E.g. it would also return 
  ## TRUE for a vector that contains a sequence of NAs at the beginning
  ## and/or end of the vector. 
{
  dists.between <- unique(diff(which(!is.na(x))))
  ##value<< logical: whether the series contains gaps or not.
  return((length(dists.between) == 1 && dists.between == 1))
}
