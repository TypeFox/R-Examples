`summary.homals` <-
function(object, ...)
{
  cat("\nNumber of dimensions:",object$ndim)
  cat("\nNumber of iterations:",object$niter)
  cat("\n")
  nvar <- dim(object$dframe)[2]
  for (i in 1:nvar) 
  {
    cat("\n---------\n")
    cat("\nVariable:",colnames(object$dframe[i]))
    cat("\nLoadings:\n")
    print(round(object$loadings[[i]], 4))
    cat("\nCategory centroids:\n")
    print(round(object$cat.centroids[[i]], 4))
    cat("\nCategory quantifications (scores):\n")
    print(round(object$catscores[[i]], 4)) 
    cat("\nLower rank quantifications (rank = ",object$rank.vec[i],"):\n",sep="")
    print(round(object$low.rank[[i]], 4)) 
  }
  cat("\n")
}

