getVecInfo <- function(
  ##title<< Compute vector summary statistics
  x ## << vector: input vector
  )
  ##description<< This function computes several summary statistics of a vector.
{
  n.na      <- sum(is.na(x)) / length(x)
  n.na.in   <- NA
  if (n.na == 1) {
    min.x     <- NA
    max.x     <- NA
    mean.x    <- NA
    sd.x      <- NA
    range.x   <- c(NA, NA)
  } else {
    range.x   <- range(x, na.rm = TRUE)
    min.x     <- range.x[1]
    max.x     <- range.x[2]
    mean.x    <- mean(x, na.rm = TRUE)
    sd.x      <- sd(x, na.rm = TRUE)
    n.na.in <- sum(is.na(x[min(which(!is.na(x))):max(which(!is.na(x)))])) / length(x)
  }
  n.inf     <- sum(is.infinite(x)) / length(x)
  ##value<< vector with the vector statistics (min,  max, mean,  sdev, range, ratio na, ratio na inner,  ratio INF)
  out       <- c(min.x, max.x, mean.x, sd.x, diff(range.x), n.na, n.na.in, n.inf)
  names(out)<- c('min', 'max', 'mean', 'sdev', 'range', 'ratio na', 'ratio na inner', 'ratio inf')
  return(out)
}  
