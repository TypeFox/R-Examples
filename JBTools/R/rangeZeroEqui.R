rangeZeroEqui <- function(
  ##title<< Compute a zero centered equi-sided range
  x ##<< input vector
  )
  ##description<< rangeZeroEqui computes a zero centered equi-sided range,  i.e. it
  ##              the maximum absolute value and returns this as a positive and a
  ##              negative value. 
  {
  ##value<< vector
  return(c(-1,1)*max(abs(range(x, na.rm = TRUE))))
}

