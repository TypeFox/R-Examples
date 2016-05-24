findCov <- function(mat,Mats,Out)
{
  # Define LISREL name:
  if (mat=="ObsCov") lisName <- "Covariance Matrix"
  if (mat=="ImpCov") lisName <- "Fitted Covariance Matrix"
  
  Res <- list()
  
  if (length(Mats[[mat]]) > 0)
  {
    for (g in seq_along(Mats[[mat]]))
    {
      Inds <- matRange(Mats[[mat]][[g]],lisName,Out)
      
      Res[[g]] <- getMatrix(Out[Inds[1]:Inds[2]],lisName,TRUE,TRUE)
    }
  }
  return(Res)
}