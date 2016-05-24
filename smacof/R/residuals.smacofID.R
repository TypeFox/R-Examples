#residuals of distances for smacof objects

residuals.smacofID <- function(object, ...)
{
  reslist <- list(NULL)
  for (i in 1:length(object$dhat)) reslist[[i]] <- (as.matrix(object$dhat[[i]] - object$confdiss[[i]]))
  names(reslist) <- names(object$dhat)
  return(reslist)  
}
