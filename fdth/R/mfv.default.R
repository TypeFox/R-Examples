mfv.default <- function(x, ...)
{
  ## Adapted from
  ## http://stackoverflow.com/questions/2547402/standard-library-function-in-r-for-finding-the-mode 
  ux <- unique(x)

  res <- as.vector(ux[which.max(tabulate(match(x, 
                                               ux)))])
  return(res)
}                        
