`segregationRatios` <-
function(x, drop.cols=NULL)
{
  
  ## Description: Computes segregation ratios for a matrix of markers
  ## where the rows are markers and the columns are individuals

  ## Arguments:
  ## x : matrix of 0's, 1,s and NA's representing scores of dominant markers
  ##     where the rows are markers and the columns are individuals
  ## drop.cols: columns to drop when calculating segregation ratios

  ## Values:
  ## Returns an object of type segRatio containing
  ## r:  no. of 1's for each individual
  ## n:  total no. of markers present for each individual
  ## seg.ratio:  segregation proportion for each individual
  ## n.individuals: total number of individuals

  if (! is.matrix(x))
    stop("A matrix must be supplied")

  r <- rowSums(x,na.rm=TRUE)
  n <- rowSums(!is.na(x))

  res <- list(r=r, n=n, seg.ratio=r/n, n.individuals=dim(x)[2],
              call=match.call())
  oldClass(res) <- "segRatio"
  return(res)
  
}

