kdtree_build <- function(data)
{
  if (class(data) == "data.frame")
  {
    data <- as.matrix(data)
  }
  
  if (class(data) == "matrix")
  {
    kdt <- kdtree_build_intl(t(data),nrow(data),ncol(data))
    return(kdt)
  }
  else
  {
    stop("Input data not a matrix or data frame")
  }
}