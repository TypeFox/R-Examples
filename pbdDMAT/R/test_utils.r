matrix.rand <- function(nrows=1, ncols=1, fill="dense", storage="dense", ..., sparsity=.1)
{
  if (storage == "dense")
  {
    if (fill == "sparse")
    {
      ind <- c(TRUE, FALSE)
      ret <- matrix(0.0, nrows, ncols)
      
      for (i in 1:ncols)
      {
        for (j in 1:nrows)
        {
          if (sample(ind, size=1, prob=c(sparsity, 1-sparsity)))
            ret[i, j] <- stats::rnorm(1)
        }
      }
    }
    else if (fill == "dense")
    {
      ret <- matrix(stats::rnorm(nrows*ncols), nrows, ncols)
    }
  }
  
  else if (storage == "sparse")
  {
    #...
  }
  
  return( ret )
}
